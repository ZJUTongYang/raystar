#include "raystar_node.h"
#include "conservative_path_length.h"
#include "metric_bound_search.h"
#include "published_path_order.h"

#include <std_msgs/msg/color_rgba.hpp>
#include <raystar_interfaces/msg/debug_node.hpp>
#include <raystar_interfaces/msg/path_result.hpp>
#include <raystar_interfaces/msg/planning_result_info.hpp>
#include <rcl_interfaces/msg/integer_range.hpp>
#include <rcl_interfaces/msg/parameter_descriptor.hpp>
#include <rcl_interfaces/msg/set_parameters_result.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <type_traits>
#include <utility>
#include <exception>

#include "raystar_node_detail.h"

namespace raystar {

using namespace node_impl;

template <typename RequestT, typename ResponseT>
void RaystarNode::executePlanning(const RequestT& request_value,
                                  ResponseT& response_value,
                                  const nav_msgs::msg::OccupancyGrid& grid,
                                  const raystar_interfaces::MapId& map_id,
                                  const StopPredicate& stop_requested) noexcept try {
  // The Core search tree and retained visualization cache are one stateful
  // unit.  The timer uses try_lock(), so it skips a tick instead of delaying
  // planning or Action cancellation while this lock is held.
  std::unique_lock<std::mutex> planner_lock(planner_cache_mutex_);
  const auto* request = &request_value;
  auto* response = &response_value;
  initializePlanningResponse(*response,
                             request->search_mode,
                             request->k,
                             request->max_path_length,
                             map_id,
                             request->include_debug);
  if (stop_requested && stop_requested()) {
    markCancelled(*response);
    response->message = "Planning was cancelled before execution started";
    return;
  }
  // Validation failures remain machine-readable even though the individual
  // branches retain precise human diagnostics below.
  response->result_info.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
  // Use only the previously admitted frame for the clear snapshot.  Copying
  // an untrusted request frame here would multiply an arbitrarily large client
  // string before the resource parameters have even been loaded.
  clearVisualizationsLocked();

  RequestConfiguration configuration;
  std::string configuration_error;
  if (!loadRequestConfiguration(configuration, configuration_error)) {
    response->success = false;
    response->result_info.status = PlanningResultInfo::STATUS_INVALID_CONFIGURATION;
    response->message = "Invalid Raystar server configuration: " + configuration_error;
    return;
  }
  response->result_info.environment_id = raystar_interfaces::computeEnvironmentId(
    grid, configuration.occupied_threshold, request->allow_unknown);
  configuration.planning.cancel_requested = stop_requested;
  if (grid.header.frame_id.size() > kMaxFrameIdBytes ||
      request->start.header.frame_id.size() > kMaxFrameIdBytes ||
      request->goal.header.frame_id.size() > kMaxFrameIdBytes) {
    response->success = false;
    response->message = "Invalid request: frame_id must be at most 256 bytes";
    return;
  }
  if (grid.header.frame_id.empty()) {
    response->success = false;
    response->message = "Invalid map: header.frame_id must not be empty";
    return;
  }
  const bool top_k_mode = request->search_mode == PlanAction::Goal::SEARCH_MODE_TOP_K;
  const bool cost_bounded_mode =
    request->search_mode == PlanAction::Goal::SEARCH_MODE_ALL_WITHIN_LENGTH;
  if (!top_k_mode && !cost_bounded_mode) {
    response->success = false;
    response->message = "Invalid search_mode: expected TOP_K or ALL_WITHIN_LENGTH";
    return;
  }
  if (top_k_mode) {
    if (request->k <= 0) {
      response->success = false;
      response->message = "Invalid K: must be greater than zero";
      return;
    }
    if (request->k > configuration.planning.max_k) {
      response->success = false;
      response->message = "Invalid K: requested " + std::to_string(request->k) +
                          " exceeds max_k=" + std::to_string(configuration.planning.max_k);
      return;
    }
    if (request->max_path_length != 0.0) {
      response->success = false;
      response->message = "Invalid TOP_K request: max_path_length must be zero";
      return;
    }
  } else {
    if (request->k != 0) {
      response->success = false;
      response->message = "Invalid ALL_WITHIN_LENGTH request: K must be zero";
      return;
    }
    if (!std::isfinite(request->max_path_length) || request->max_path_length <= 0.0) {
      response->success = false;
      response->message =
        "Invalid ALL_WITHIN_LENGTH request: max_path_length must be finite and greater than zero";
      return;
    }
  }

  // Validate the small, scalar endpoint fields before copying the map.  This
  // keeps an otherwise rejected request from paying the binary-map allocation
  // cost (the frame-length checks above provide the same guard for strings).
  std::string pose_error;
  if (!validatePlanarPose(request->start, "Start", grid.header.frame_id, pose_error)) {
    response->success = false;
    response->message = pose_error;
    return;
  }
  if (!validatePlanarPose(request->goal, "Goal", grid.header.frame_id, pose_error)) {
    response->success = false;
    response->message = pose_error;
    return;
  }

  GridMap work_map;
  std::string map_error;
  if (!occupancyGridToBinaryMap(
        grid, request->allow_unknown, configuration, stop_requested, work_map, map_error)) {
    response->success = false;
    if (stop_requested && stop_requested())
      markCancelled(*response);
    response->message = map_error;
    RCLCPP_WARN(get_logger(), "Rejected OccupancyGrid: %s", map_error.c_str());
    return;
  }

  ContinuousGridPoint start_point;
  ContinuousGridPoint goal_point;
  if (!worldToContinuousMap(
        work_map, request->start.pose.position.x, request->start.pose.position.y, start_point)) {
    response->success = false;
    response->message = "Start position is outside the map";
    return;
  }
  if (!worldToContinuousMap(
        work_map, request->goal.pose.position.x, request->goal.pose.position.y, goal_point)) {
    response->success = false;
    response->message = "Goal position is outside the map";
    return;
  }

  const PlanEndpoint start_endpoint(
    {static_cast<int>(start_point.cell_x), static_cast<int>(start_point.cell_y)},
    {start_point.x, start_point.y});
  const PlanEndpoint goal_endpoint(
    {static_cast<int>(goal_point.cell_x), static_cast<int>(goal_point.cell_y)},
    {goal_point.x, goal_point.y});

  SearchObjective objective = SearchObjective::topK(request->k);
  if (cost_bounded_mode) {
    double padded_metric_bound = std::numeric_limits<double>::quiet_NaN();
    if (!paddedMetricBoundForGridSearch(work_map,
                                        request->max_path_length,
                                        configuration.planning.max_nodes,
                                        padded_metric_bound)) {
      response->success = false;
      response->message =
        "Invalid max_path_length: could not bound grid-to-world representation error";
      return;
    }
    double grid_cost = std::numeric_limits<double>::quiet_NaN();
    if (!gridCostSearchUpperBoundForMetricBound(
          padded_metric_bound, work_map.resolution, grid_cost)) {
      response->success = false;
      response->message =
        "Invalid max_path_length: could not derive its binary64 Core search bound";
      return;
    }
    const auto metric_admission = [&work_map, metric_bound = request->max_path_length](
                                    const BoundedPathView& path, const StopToken& stop_token) {
      switch (classifyPathViewMetricBound(path, work_map, metric_bound, stop_token)) {
      case MetricBoundEligibility::within_bound:
        return BoundedPathAdmission::accept;
      case MetricBoundEligibility::outside_bound:
        return BoundedPathAdmission::reject;
      case MetricBoundEligibility::invalid:
        return BoundedPathAdmission::failure_or_stop;
      }
      return BoundedPathAdmission::failure_or_stop;
    };
    objective = SearchObjective::allWithinCost(grid_cost, metric_admission);
  }

  RCLCPP_INFO(get_logger(),
              "Planning: start=(%.6f,%.6f) cell=(%u,%u) goal=(%.6f,%.6f) "
              "cell=(%u,%u) mode=%s K=%d max_path_length=%.9f "
              "allow_self_crossing=%s allow_unknown=%s "
              "occupied_threshold=%d",
              start_point.x,
              start_point.y,
              start_point.cell_x,
              start_point.cell_y,
              goal_point.x,
              goal_point.y,
              goal_point.cell_x,
              goal_point.cell_y,
              top_k_mode ? "TOP_K" : "ALL_WITHIN_LENGTH",
              request->k,
              request->max_path_length,
              request->allow_self_crossing ? "true" : "false",
              request->allow_unknown ? "true" : "false",
              configuration.occupied_threshold);

  auto result = core_.plan(work_map,
                           start_endpoint,
                           goal_endpoint,
                           objective,
                           request->allow_self_crossing,
                           configuration.planning);
  SearchStateRelease search_state_release(core_);
  const auto& nodes = core_.getNodes();

  auto& result_info = response->result_info;
  result_info.found_path_count = boundedUint32(result.path_solutions.size());
  result_info.expanded_nodes = boundedUint64(result.expanded_nodes);
  result_info.map_time_ms = result.map_time_ms;
  result_info.plan_time_ms = result.plan_time_ms;
  result_info.limits_reached = planningLimitMask(result.limit_reached);
  result_info.cost_bound_exhausted =
    cost_bounded_mode && (result.completion == PlanningCompletion::cost_bound_exhausted ||
                          result.completion == PlanningCompletion::frontier_exhausted);
  switch (result.outcome) {
  case PlanningOutcome::complete:
    result_info.search_complete = true;
    result_info.status =
      cost_bounded_mode || result.path_solutions.size() >= static_cast<size_t>(request->k)
        ? PlanningResultInfo::STATUS_COMPLETE
        : PlanningResultInfo::STATUS_FEWER_PATHS;
    break;
  case PlanningOutcome::no_path:
    result_info.search_complete = true;
    result_info.status = PlanningResultInfo::STATUS_NO_PATH;
    break;
  case PlanningOutcome::limit_reached:
    result_info.search_complete = false;
    result_info.status = result.limit_reached == PlanningLimitReached::cancelled
                           ? PlanningResultInfo::STATUS_CANCELLED
                           : PlanningResultInfo::STATUS_PARTIAL_SEARCH;
    break;
  case PlanningOutcome::invalid_request:
    result_info.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
    break;
  case PlanningOutcome::failed:
    result_info.status = PlanningResultInfo::STATUS_FAILED;
    break;
  }

  RCLCPP_INFO(get_logger(),
              "Planning done: map=%.1fms plan=%.1fms paths=%zu",
              result.map_time_ms,
              result.plan_time_ms,
              result.path_solutions.size());

  if (stop_requested && stop_requested()) {
    markCancelled(*response);
    response->message = result.message.empty() ? "Planning was cancelled" : result.message;
    return;
  }

  response->message = result.message;
  if (response->message.size() > kMaxDiagnosticBytes)
    response->message.resize(kMaxDiagnosticBytes);

  const double resolution = static_cast<double>(work_map.resolution);
  size_t response_bytes = kResponseBaseBytes;
  if (!checkedAdd(
        response_bytes, response->message.size(), configuration.planning.max_response_bytes) ||
      !checkedAdd(
        response_bytes, grid.header.frame_id.size(), configuration.planning.max_response_bytes)) {
    response->success = false;
    markPathOutputIncomplete(result_info, PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
    result_info.returned_path_count = 0;
    result_info.request_satisfied = false;
    response->message = "Response metadata exceeds max_response_bytes=" +
                        std::to_string(configuration.planning.max_response_bytes);
    return;
  }

  // Keep only indices while serializing.  Copying every accepted PathSolution
  // here would duplicate turning-point/node vectors on top of Core's result;
  // the selected solutions are moved into the visualization cache only after
  // the response has passed its output limits.
  std::vector<size_t> metric_eligible_solution_indices;
  metric_eligible_solution_indices.reserve(result.path_solutions.size());
  size_t precertification_omitted_count = 0;
  for (size_t solution_index = 0; solution_index < result.path_solutions.size(); ++solution_index) {
    if (!cost_bounded_mode) {
      metric_eligible_solution_indices.emplace_back(solution_index);
      continue;
    }
    if (stop_requested && stop_requested()) {
      initializePlanningResponse(*response,
                                 request->search_mode,
                                 request->k,
                                 request->max_path_length,
                                 map_id,
                                 request->include_debug);
      result_info.found_path_count = boundedUint32(result.path_solutions.size());
      result_info.expanded_nodes = boundedUint64(result.expanded_nodes);
      result_info.map_time_ms = result.map_time_ms;
      result_info.plan_time_ms = result.plan_time_ms;
      markCancelled(*response);
      response->message = "Planning was cancelled while certifying topology output";
      return;
    }
    nav_msgs::msg::Path topology_path;
    std::string topology_error;
    const auto& solution = result.path_solutions[solution_index];
    if (!buildTopologyPathMsg(
          solution, work_map, grid.header.frame_id, topology_path, topology_error)) {
      ++precertification_omitted_count;
      appendResponseNotice(response->message,
                           "A topology path was omitted: " + topology_error,
                           response_bytes,
                           configuration.planning.max_response_bytes);
      continue;
    }
    const double core_metric_cost = gridCostToMetric(solution.path_cost_, work_map.resolution);
    const auto eligibility = classifyTopologyMetricBound(
      topology_path, core_metric_cost, request->max_path_length, stop_requested);
    if (eligibility == MetricBoundEligibility::outside_bound) {
      continue;
    }
    if (eligibility == MetricBoundEligibility::invalid) {
      ++precertification_omitted_count;
      appendResponseNotice(
        response->message,
        "A topology path was omitted because its metric length could not be certified",
        response_bytes,
        configuration.planning.max_response_bytes);
      continue;
    }
    metric_eligible_solution_indices.emplace_back(solution_index);
  }
  result_info.found_path_count = boundedUint32(metric_eligible_solution_indices.size());

  std::vector<size_t> emitted_solution_indices;
  size_t omitted_path_count = precertification_omitted_count;
  size_t emitted_path_points = 0;
  uint16_t path_output_limits = PlanningResultInfo::LIMIT_NONE;
  for (const size_t solution_index : metric_eligible_solution_indices) {
    if (stop_requested && stop_requested()) {
      initializePlanningResponse(*response,
                                 request->search_mode,
                                 request->k,
                                 request->max_path_length,
                                 map_id,
                                 request->include_debug);
      result_info.found_path_count = boundedUint32(result.path_solutions.size());
      result_info.expanded_nodes = boundedUint64(result.expanded_nodes);
      result_info.map_time_ms = result.map_time_ms;
      result_info.plan_time_ms = result.plan_time_ms;
      markCancelled(*response);
      response->message = "Planning was cancelled while building the response";
      return;
    }
    const auto& sol = result.path_solutions[solution_index];
    raystar_interfaces::msg::PathResult path_result;
    std::string path_error;
    if (!buildTopologyPathMsg(
          sol, work_map, grid.header.frame_id, path_result.topology_path, path_error)) {
      ++omitted_path_count;
      appendResponseNotice(response->message,
                           "A topology path was omitted: " + path_error,
                           response_bytes,
                           configuration.planning.max_response_bytes);
      continue;
    }
    const size_t topology_point_count = path_result.topology_path.poses.size();
    size_t point_count = 0;
    bool topology_only = false;
    if (!countInterpolatedPathPoints(
          sol, configuration.planning.max_path_points, point_count, path_error)) {
      if (cost_bounded_mode) {
        topology_only = true;
        point_count = topology_point_count;
      } else {
        ++omitted_path_count;
        path_output_limits =
          static_cast<uint16_t>(path_output_limits | PlanningResultInfo::LIMIT_MAX_PATH_POINTS);
        appendResponseNotice(response->message,
                             "A path was omitted: " + path_error,
                             response_bytes,
                             configuration.planning.max_response_bytes);
        continue;
      }
    }
    const auto points_fit = [&](size_t serialized_path_points) {
      return topology_point_count <= configuration.planning.max_path_points - emitted_path_points &&
             serialized_path_points <=
               configuration.planning.max_path_points - emitted_path_points - topology_point_count;
    };
    if (!points_fit(point_count) && cost_bounded_mode && !topology_only &&
        points_fit(topology_point_count)) {
      topology_only = true;
      point_count = topology_point_count;
    }
    if (!points_fit(point_count)) {
      omitted_path_count = precertification_omitted_count +
                           metric_eligible_solution_indices.size() -
                           emitted_solution_indices.size();
      path_output_limits =
        static_cast<uint16_t>(path_output_limits | PlanningResultInfo::LIMIT_MAX_PATH_POINTS);
      appendResponseNotice(response->message,
                           "Path output reached the per-response max_path_points=" +
                             std::to_string(configuration.planning.max_path_points),
                           response_bytes,
                           configuration.planning.max_response_bytes);
      break;
    }

    size_t path_bytes = 0;
    size_t topology_path_bytes = 0;
    size_t candidate_response_bytes = response_bytes;
    const auto bytes_fit = [&](size_t serialized_path_points, size_t& candidate_bytes) {
      candidate_bytes = response_bytes;
      return estimatePathResponseBytes(
               serialized_path_points, grid.header.frame_id.size(), path_bytes) &&
             estimatePathResponseBytes(
               topology_point_count, grid.header.frame_id.size(), topology_path_bytes) &&
             checkedAdd(candidate_bytes, path_bytes, configuration.planning.max_response_bytes) &&
             checkedAdd(
               candidate_bytes, topology_path_bytes, configuration.planning.max_response_bytes);
    };
    bool response_bytes_fit = bytes_fit(point_count, candidate_response_bytes);
    if (!response_bytes_fit && cost_bounded_mode && !topology_only) {
      topology_only = true;
      point_count = topology_point_count;
      response_bytes_fit = bytes_fit(point_count, candidate_response_bytes);
    }
    if (!response_bytes_fit) {
      path_output_limits =
        static_cast<uint16_t>(path_output_limits | PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
      appendResponseNotice(response->message,
                           "Path output reached max_response_bytes=" +
                             std::to_string(configuration.planning.max_response_bytes),
                           response_bytes,
                           configuration.planning.max_response_bytes);
      // Paths are ordered by Core's cost/topology ranking.  Once one no longer
      // fits, later paths cannot make the response smaller in a predictable
      // way, so stop before allocating another Path message.
      omitted_path_count = precertification_omitted_count +
                           metric_eligible_solution_indices.size() -
                           emitted_solution_indices.size();
      break;
    }

    if (topology_only) {
      path_result.path = path_result.topology_path;
    } else if (!buildPathMsg(sol,
                             work_map,
                             grid.header.frame_id,
                             configuration.planning.max_path_points - emitted_path_points -
                               topology_point_count,
                             path_result.path,
                             path_error)) {
      ++omitted_path_count;
      appendResponseNotice(response->message,
                           "A path was omitted: " + path_error,
                           response_bytes,
                           configuration.planning.max_response_bytes);
      continue;
    }
    const double core_metric_cost = gridCostToMetric(sol.path_cost_, work_map.resolution);
    const bool finalized =
      cost_bounded_mode
        ? finalizeCostBoundedPublishedPathResult(
            path_result, core_metric_cost, request->max_path_length, stop_requested)
        : finalizePublishedPathResult(
            path_result, core_metric_cost, std::numeric_limits<double>::infinity(), stop_requested);
    if (!finalized) {
      ++omitted_path_count;
      appendResponseNotice(response->message,
                           "A path was omitted because its serialized world geometry could not "
                           "faithfully represent the bounded Core cost",
                           response_bytes,
                           configuration.planning.max_response_bytes);
      continue;
    }
    // The bounded finalizer may replace a dense interpolation with the
    // shorter topology path at a one-ULP metric boundary. Charge the geometry
    // actually serialized, not the earlier dense estimate, so one fallback
    // cannot spuriously consume another path's response budget.
    if (!bytes_fit(path_result.path.poses.size(), candidate_response_bytes)) {
      ++omitted_path_count;
      path_output_limits =
        static_cast<uint16_t>(path_output_limits | PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
      appendResponseNotice(response->message,
                           "A certified path did not fit max_response_bytes=" +
                             std::to_string(configuration.planning.max_response_bytes),
                           response_bytes,
                           configuration.planning.max_response_bytes);
      continue;
    }
    const size_t serialized_path_points =
      path_result.path.poses.size() + path_result.topology_path.poses.size();
    response_bytes = candidate_response_bytes;
    response->path_results.emplace_back(std::move(path_result));
    emitted_solution_indices.push_back(solution_index);
    emitted_path_points += serialized_path_points;
  }

  if (!stableSortPublishedPathsWithSources(response->path_results, emitted_solution_indices)) {
    throw std::logic_error("published path/source index cardinality mismatch");
  }

  if (omitted_path_count > 0) {
    appendResponseNotice(
      response->message,
      std::to_string(omitted_path_count) + " path(s) omitted while building ROS output",
      response_bytes,
      configuration.planning.max_response_bytes);
  }
  result_info.returned_path_count = boundedUint32(response->path_results.size());
  const size_t metric_found_path_count = metric_eligible_solution_indices.size();
  result_info.output_complete =
    omitted_path_count == 0 && response->path_results.size() == metric_found_path_count;
  if (!result_info.output_complete)
    markPathOutputIncomplete(result_info, path_output_limits);
  result_info.request_satisfied =
    cost_bounded_mode ? result_info.cost_bound_exhausted && result_info.output_complete
                      : result_info.requested_path_count > 0 &&
                          result_info.returned_path_count == result_info.requested_path_count;
  response->success = !response->path_results.empty();
  if (cost_bounded_mode && result_info.search_complete && metric_found_path_count == 0)
    result_info.status = PlanningResultInfo::STATUS_NO_PATH;

  const size_t response_budget_remaining =
    configuration.planning.max_response_bytes > response_bytes
      ? configuration.planning.max_response_bytes - response_bytes
      : 0;
  const size_t debug_by_bytes = response_budget_remaining / kDebugNodeBytes;
  const size_t debug_node_count =
    request->include_debug
      ? std::min({nodes.size(), configuration.planning.max_debug_nodes, debug_by_bytes})
      : 0U;
  response->debug_nodes.reserve(debug_node_count);
  for (size_t i = 0; i < debug_node_count; ++i) {
    if ((i & 0xffu) == 0u && stop_requested && stop_requested()) {
      initializePlanningResponse(*response,
                                 request->search_mode,
                                 request->k,
                                 request->max_path_length,
                                 map_id,
                                 request->include_debug);
      result_info.found_path_count = boundedUint32(result.path_solutions.size());
      result_info.expanded_nodes = boundedUint64(result.expanded_nodes);
      result_info.map_time_ms = result.map_time_ms;
      result_info.plan_time_ms = result.plan_time_ms;
      markCancelled(*response);
      response->message = "Planning was cancelled while building debug output";
      return;
    }
    const auto& n = nodes[i];
    raystar_interfaces::msg::DebugNode debug_node;
    debug_node.index = n.index();
    debug_node.g_cost = n.gCost() * resolution;
    debug_node.f_cost = (n.gCost() + n.hCost()) * resolution;
    response->debug_nodes.emplace_back(std::move(debug_node));
  }
  const size_t emitted_debug_nodes = debug_node_count;
  // debug_node_count was admitted from response_budget_remaining /
  // kDebugNodeBytes, so both the multiplication and addition are bounded.
  // Charge it before appending a truncation notice; otherwise the notice can
  // consume bytes already occupied by serialized debug records.
  response_bytes += emitted_debug_nodes * kDebugNodeBytes;
  if (request->include_debug && emitted_debug_nodes < result.expanded_nodes) {
    appendResponseNotice(
      response->message,
      "Debug output limited to " + std::to_string(emitted_debug_nodes) + " node(s)",
      response_bytes,
      configuration.planning.max_response_bytes);
  }

  result_info.debug_output_complete =
    !request->include_debug || emitted_debug_nodes == result.expanded_nodes;
  response->success = !response->path_results.empty();

  if (stop_requested && stop_requested()) {
    const auto environment_id = response->result_info.environment_id;
    initializePlanningResponse(*response,
                               request->search_mode,
                               request->k,
                               request->max_path_length,
                               map_id,
                               request->include_debug);
    response->result_info.environment_id = environment_id;
    result_info.found_path_count = boundedUint32(result.path_solutions.size());
    result_info.expanded_nodes = boundedUint64(result.expanded_nodes);
    result_info.map_time_ms = result.map_time_ms;
    result_info.plan_time_ms = result.plan_time_ms;
    markCancelled(*response);
    response->message = "Planning was cancelled before publishing results";
    return;
  }

  if (response->success && result.polymap) {
    std::vector<PathSolution> visualization_solutions;
    visualization_solutions.reserve(emitted_solution_indices.size());
    for (const size_t solution_index : emitted_solution_indices) {
      visualization_solutions.emplace_back(std::move(result.path_solutions[solution_index]));
    }
    last_frame_id_ = grid.header.frame_id;

    publishNonHomotopicPaths(visualization_solutions,
                             work_map,
                             grid.header.frame_id,
                             configuration.planning.max_path_points,
                             configuration.planning.max_response_bytes);
    if (request->include_debug) {
      publishPolyObstacles(
        *result.polymap, work_map, grid.header.frame_id, configuration.planning.max_response_bytes);
      publishCDT(
        *result.polymap, work_map, grid.header.frame_id, configuration.planning.max_response_bytes);
      publishDebugTree(nodes,
                       work_map,
                       grid.header.frame_id,
                       configuration.planning.max_debug_nodes,
                       configuration.planning.max_response_bytes);
    }
  }

} catch (const std::exception& exception) {
  std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
  clearVisualizationsLocked();
  initializePlanningResponse(response_value,
                             request_value.search_mode,
                             request_value.k,
                             request_value.max_path_length,
                             map_id,
                             request_value.include_debug);
  if (stop_requested && stop_requested()) {
    markCancelled(response_value);
    response_value.message = "Planning was cancelled while handling an exception";
  } else {
    response_value.result_info.status = PlanningResultInfo::STATUS_FAILED;
    setBoundedExceptionMessage(response_value.message, exception.what());
  }
  RCLCPP_ERROR(get_logger(), "%s", response_value.message.c_str());
} catch (...) {
  std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
  clearVisualizationsLocked();
  initializePlanningResponse(response_value,
                             request_value.search_mode,
                             request_value.k,
                             request_value.max_path_length,
                             map_id,
                             request_value.include_debug);
  if (stop_requested && stop_requested()) {
    markCancelled(response_value);
    response_value.message = "Planning was cancelled while handling an unknown exception";
  } else {
    response_value.result_info.status = PlanningResultInfo::STATUS_FAILED;
    response_value.message = "Raystar request failed with an unknown exception";
  }
  RCLCPP_ERROR(get_logger(), "%s", response_value.message.c_str());
}

void RaystarNode::executeGoalSetPlanning(const GoalSetAction::Goal& request,
                                         GoalSetAction::Result& response,
                                         const nav_msgs::msg::OccupancyGrid& grid,
                                         const raystar_interfaces::MapId& map_id,
                                         const StopPredicate& stop_requested) noexcept try {
  std::unique_lock<std::mutex> planner_lock(planner_cache_mutex_);
  response = GoalSetAction::Result{};
  response.requested_start = request.start;
  response.requested_allow_self_crossing = request.allow_self_crossing;
  response.requested_allow_unknown = request.allow_unknown;
  response.result_info.map_id = map_id;
  response.result_info.requested_goal_count = boundedUint32(request.goals.size());
  response.result_info.debug_requested = request.include_debug;
  response.result_info.debug_output_complete = !request.include_debug;
  response.result_info.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
  clearVisualizationsLocked();
  const auto mark_goal_set_cancelled = [&]() {
    response.success = false;
    response.goal_results.clear();
    response.debug_nodes.clear();
    response.result_info.status = PlanningResultInfo::STATUS_CANCELLED;
    response.result_info.limits_reached = static_cast<uint16_t>(
      response.result_info.limits_reached | PlanningResultInfo::LIMIT_CANCELLED);
    response.result_info.request_satisfied = false;
    response.result_info.search_complete = false;
    response.result_info.output_complete = false;
    response.result_info.debug_output_complete = false;
    response.result_info.returned_goal_count = 0;
    response.result_info.completed_goal_count = 0;
    response.result_info.goals_with_paths = 0;
    response.result_info.found_path_count = 0;
    response.result_info.returned_path_count = 0;
    response.message = "Planning was cancelled while certifying goal-set output";
  };

  RequestConfiguration configuration;
  std::string configuration_error;
  if (!loadRequestConfiguration(configuration, configuration_error)) {
    response.result_info.status = PlanningResultInfo::STATUS_INVALID_CONFIGURATION;
    response.message = "Invalid Raystar server configuration: " + configuration_error;
    return;
  }
  response.result_info.environment_id = raystar_interfaces::computeEnvironmentId(
    grid, configuration.occupied_threshold, request.allow_unknown);
  configuration.planning.cancel_requested = stop_requested;
  if (request.goals.empty()) {
    response.message = "Invalid goal-set request: at least one goal is required";
    return;
  }
  if (request.goals.size() != request.max_path_lengths.size()) {
    response.message =
      "Invalid goal-set request: goals and max_path_lengths must have equal lengths";
    return;
  }
  if (request.goals.size() > configuration.planning.max_multi_goal_count) {
    response.message = "Invalid goal-set request: requested " +
                       std::to_string(request.goals.size()) +
                       " goals exceeds max_multi_goal_count=" +
                       std::to_string(configuration.planning.max_multi_goal_count);
    return;
  }
  if (grid.header.frame_id.empty()) {
    response.message = "Invalid map: header.frame_id must not be empty";
    return;
  }
  if (grid.header.frame_id.size() > kMaxFrameIdBytes ||
      request.start.header.frame_id.size() > kMaxFrameIdBytes) {
    response.message = "Invalid request: frame_id must be at most 256 bytes";
    return;
  }
  for (size_t i = 0; i < request.goals.size(); ++i) {
    if (request.goals[i].header.frame_id.size() > kMaxFrameIdBytes) {
      response.message =
        "Invalid goal[" + std::to_string(i) + "]: frame_id must be at most 256 bytes";
      return;
    }
    if (!std::isfinite(request.max_path_lengths[i]) || request.max_path_lengths[i] <= 0.0) {
      response.message = "Invalid max_path_lengths[" + std::to_string(i) +
                         "]: value must be finite and greater than zero";
      return;
    }
  }

  std::string pose_error;
  if (!validatePlanarPose(request.start, "Start", grid.header.frame_id, pose_error)) {
    response.message = pose_error;
    return;
  }
  for (size_t i = 0; i < request.goals.size(); ++i) {
    const std::string label = "Goal[" + std::to_string(i) + "]";
    if (!validatePlanarPose(request.goals[i], label.c_str(), grid.header.frame_id, pose_error)) {
      response.message = pose_error;
      return;
    }
  }

  GridMap work_map;
  std::string map_error;
  if (!occupancyGridToBinaryMap(
        grid, request.allow_unknown, configuration, stop_requested, work_map, map_error)) {
    response.message = map_error;
    response.result_info.status = stop_requested && stop_requested()
                                    ? PlanningResultInfo::STATUS_CANCELLED
                                    : PlanningResultInfo::STATUS_INVALID_REQUEST;
    return;
  }

  ContinuousGridPoint start_point;
  if (!worldToContinuousMap(
        work_map, request.start.pose.position.x, request.start.pose.position.y, start_point)) {
    response.message = "Start position is outside the map";
    return;
  }
  const PlanEndpoint start_endpoint(
    {static_cast<int>(start_point.cell_x), static_cast<int>(start_point.cell_y)},
    {start_point.x, start_point.y});
  std::vector<CostBoundedGoal> core_goals;
  core_goals.reserve(request.goals.size());
  for (size_t i = 0; i < request.goals.size(); ++i) {
    ContinuousGridPoint goal_point;
    if (!worldToContinuousMap(work_map,
                              request.goals[i].pose.position.x,
                              request.goals[i].pose.position.y,
                              goal_point)) {
      response.message = "Goal[" + std::to_string(i) + "] position is outside the map";
      return;
    }
    double padded_metric_bound = std::numeric_limits<double>::quiet_NaN();
    if (!paddedMetricBoundForGridSearch(work_map,
                                        request.max_path_lengths[i],
                                        configuration.planning.max_nodes,
                                        padded_metric_bound)) {
      response.message = "Invalid max_path_lengths[" + std::to_string(i) +
                         "]: could not bound grid-to-world representation error";
      return;
    }
    double grid_cost = std::numeric_limits<double>::quiet_NaN();
    if (!gridCostSearchUpperBoundForMetricBound(
          padded_metric_bound, work_map.resolution, grid_cost)) {
      response.message = "Invalid max_path_lengths[" + std::to_string(i) +
                         "]: could not derive its binary64 Core search bound";
      return;
    }
    const auto metric_admission = [&work_map, metric_bound = request.max_path_lengths[i]](
                                    const BoundedPathView& path, const StopToken& stop_token) {
      switch (classifyPathViewMetricBound(path, work_map, metric_bound, stop_token)) {
      case MetricBoundEligibility::within_bound:
        return BoundedPathAdmission::accept;
      case MetricBoundEligibility::outside_bound:
        return BoundedPathAdmission::reject;
      case MetricBoundEligibility::invalid:
        return BoundedPathAdmission::failure_or_stop;
      }
      return BoundedPathAdmission::failure_or_stop;
    };
    core_goals.push_back(
      {PlanEndpoint({static_cast<int>(goal_point.cell_x), static_cast<int>(goal_point.cell_y)},
                    {goal_point.x, goal_point.y}),
       grid_cost,
       metric_admission});
  }

  RCLCPP_INFO(get_logger(),
              "Multi-goal planning: start=(%.6f,%.6f) goals=%zu allow_self_crossing=%s "
              "allow_unknown=%s occupied_threshold=%d",
              start_point.x,
              start_point.y,
              core_goals.size(),
              request.allow_self_crossing ? "true" : "false",
              request.allow_unknown ? "true" : "false",
              configuration.occupied_threshold);
  auto result = core_.planToGoalsWithinCosts(
    work_map, start_endpoint, core_goals, request.allow_self_crossing, configuration.planning);
  SearchStateRelease search_state_release(core_);
  const auto& nodes = core_.getNodes();
  if (stop_requested && stop_requested()) {
    mark_goal_set_cancelled();
    return;
  }

  auto& aggregate = response.result_info;
  aggregate.expanded_nodes = boundedUint64(result.expanded_nodes);
  aggregate.map_time_ms = result.map_time_ms;
  aggregate.plan_time_ms = result.plan_time_ms;
  aggregate.limits_reached = planningLimitMask(result.limit_reached);
  switch (result.outcome) {
  case PlanningOutcome::complete:
    aggregate.status = PlanningResultInfo::STATUS_COMPLETE;
    break;
  case PlanningOutcome::no_path:
    aggregate.status = PlanningResultInfo::STATUS_NO_PATH;
    break;
  case PlanningOutcome::limit_reached:
    aggregate.status = result.limit_reached == PlanningLimitReached::cancelled
                         ? PlanningResultInfo::STATUS_CANCELLED
                         : PlanningResultInfo::STATUS_PARTIAL_SEARCH;
    break;
  case PlanningOutcome::invalid_request:
    aggregate.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
    break;
  case PlanningOutcome::failed:
    aggregate.status = PlanningResultInfo::STATUS_FAILED;
    break;
  }
  response.message = result.message.substr(0, kMaxDiagnosticBytes);
  response.goal_results.reserve(result.goal_results.size());
  size_t response_bytes = kResponseBaseBytes;
  if (!checkedAdd(response_bytes,
                  response.requested_start.header.frame_id.size(),
                  configuration.planning.max_response_bytes)) {
    response.message = "Request-echo metadata exceeds max_response_bytes=" +
                       std::to_string(configuration.planning.max_response_bytes);
    aggregate.status = PlanningResultInfo::STATUS_PARTIAL_OUTPUT;
    aggregate.limits_reached = static_cast<uint16_t>(aggregate.limits_reached |
                                                     PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
    aggregate.output_complete = false;
    return;
  }
  if (!checkedAdd(
        response_bytes, response.message.size(), configuration.planning.max_response_bytes)) {
    response.message = "Response metadata exceeds max_response_bytes=" +
                       std::to_string(configuration.planning.max_response_bytes);
    aggregate.status = PlanningResultInfo::STATUS_PARTIAL_OUTPUT;
    aggregate.limits_reached = static_cast<uint16_t>(aggregate.limits_reached |
                                                     PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
    aggregate.output_complete = false;
    return;
  }
  size_t emitted_path_points = 0;
  size_t total_found_paths = 0;
  size_t total_returned_paths = 0;
  size_t completed_goals = 0;
  size_t goals_with_paths = 0;
  std::vector<std::pair<size_t, size_t>> emitted_solutions;
  bool all_output_complete = true;
  bool all_requests_satisfied = true;
  bool all_search_complete = true;
  const double resolution = static_cast<double>(work_map.resolution);

  for (const auto& core_goal_result : result.goal_results) {
    total_found_paths += core_goal_result.path_solutions.size();
    aggregate.limits_reached = static_cast<uint16_t>(
      aggregate.limits_reached | planningLimitMask(core_goal_result.limit_reached));
    const bool search_complete = core_goal_result.outcome == PlanningOutcome::complete ||
                                 core_goal_result.outcome == PlanningOutcome::no_path;
    if (search_complete)
      ++completed_goals;
    all_search_complete = all_search_complete && search_complete;
  }

  for (size_t i = 0; i < result.goal_results.size(); ++i) {
    const auto& core_goal_result = result.goal_results[i];
    raystar_interfaces::msg::GoalPathResult goal_output;
    goal_output.goal_index = boundedUint32(i);
    goal_output.goal = request.goals[i];
    goal_output.requested_max_path_length = request.max_path_lengths[i];
    goal_output.message = core_goal_result.message.substr(0, kMaxDiagnosticBytes);
    size_t goal_metadata_bytes = kGoalResultBaseBytes;
    if (!checkedAdd(goal_metadata_bytes,
                    goal_output.goal.header.frame_id.size(),
                    configuration.planning.max_response_bytes) ||
        !checkedAdd(goal_metadata_bytes,
                    goal_output.message.size(),
                    configuration.planning.max_response_bytes) ||
        !checkedAdd(
          response_bytes, goal_metadata_bytes, configuration.planning.max_response_bytes)) {
      all_output_complete = false;
      all_requests_satisfied = false;
      aggregate.limits_reached = static_cast<uint16_t>(
        aggregate.limits_reached | PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
      break;
    }
    auto& info = goal_output.result_info;
    info.map_id = map_id;
    info.environment_id = aggregate.environment_id;
    info.search_mode = PlanAction::Goal::SEARCH_MODE_ALL_WITHIN_LENGTH;
    info.requested_max_path_length = request.max_path_lengths[i];
    info.debug_requested = request.include_debug;
    info.debug_output_complete = !request.include_debug;
    info.found_path_count = boundedUint32(core_goal_result.path_solutions.size());
    info.expanded_nodes = boundedUint64(result.expanded_nodes);
    info.map_time_ms = result.map_time_ms;
    info.plan_time_ms = result.plan_time_ms;
    info.limits_reached = planningLimitMask(core_goal_result.limit_reached);
    info.cost_bound_exhausted =
      core_goal_result.completion == PlanningCompletion::cost_bound_exhausted ||
      core_goal_result.completion == PlanningCompletion::frontier_exhausted;
    switch (core_goal_result.outcome) {
    case PlanningOutcome::complete:
      info.status = PlanningResultInfo::STATUS_COMPLETE;
      info.search_complete = true;
      break;
    case PlanningOutcome::no_path:
      info.status = PlanningResultInfo::STATUS_NO_PATH;
      info.search_complete = true;
      break;
    case PlanningOutcome::limit_reached:
      info.status = core_goal_result.limit_reached == PlanningLimitReached::cancelled
                      ? PlanningResultInfo::STATUS_CANCELLED
                      : PlanningResultInfo::STATUS_PARTIAL_SEARCH;
      break;
    case PlanningOutcome::invalid_request:
      info.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
      break;
    case PlanningOutcome::failed:
      info.status = PlanningResultInfo::STATUS_FAILED;
      break;
    }

    std::vector<size_t> metric_eligible_solution_indices;
    metric_eligible_solution_indices.reserve(core_goal_result.path_solutions.size());
    size_t omitted = 0;
    for (size_t solution_index = 0; solution_index < core_goal_result.path_solutions.size();
         ++solution_index) {
      if (stop_requested && stop_requested()) {
        mark_goal_set_cancelled();
        return;
      }
      const auto& solution = core_goal_result.path_solutions[solution_index];
      nav_msgs::msg::Path topology_path;
      std::string topology_error;
      if (!buildTopologyPathMsg(
            solution, work_map, grid.header.frame_id, topology_path, topology_error)) {
        ++omitted;
        continue;
      }
      const double core_metric_cost = gridCostToMetric(solution.path_cost_, work_map.resolution);
      const auto eligibility = classifyTopologyMetricBound(
        topology_path, core_metric_cost, request.max_path_lengths[i], stop_requested);
      if (eligibility == MetricBoundEligibility::outside_bound) {
        continue;
      }
      if (eligibility == MetricBoundEligibility::invalid) {
        ++omitted;
        continue;
      }
      metric_eligible_solution_indices.emplace_back(solution_index);
    }
    const size_t metric_found_path_count = metric_eligible_solution_indices.size();
    info.found_path_count = boundedUint32(metric_found_path_count);
    std::vector<size_t> goal_emitted_solution_indices;
    goal_emitted_solution_indices.reserve(metric_found_path_count);

    for (const size_t solution_index : metric_eligible_solution_indices) {
      if (stop_requested && stop_requested()) {
        mark_goal_set_cancelled();
        return;
      }
      const auto& solution = core_goal_result.path_solutions[solution_index];
      raystar_interfaces::msg::PathResult path_result;
      std::string path_error;
      if (!buildTopologyPathMsg(
            solution, work_map, grid.header.frame_id, path_result.topology_path, path_error)) {
        ++omitted;
        continue;
      }
      const size_t topology_point_count = path_result.topology_path.poses.size();
      size_t point_count = 0;
      bool topology_only = false;
      if (!countInterpolatedPathPoints(
            solution, configuration.planning.max_path_points, point_count, path_error)) {
        topology_only = true;
        point_count = topology_point_count;
      }
      const auto points_fit = [&](size_t serialized_path_points) {
        return topology_point_count <=
                 configuration.planning.max_path_points - emitted_path_points &&
               serialized_path_points <= configuration.planning.max_path_points -
                                           emitted_path_points - topology_point_count;
      };
      if (!points_fit(point_count) && !topology_only && points_fit(topology_point_count)) {
        topology_only = true;
        point_count = topology_point_count;
      }
      if (!points_fit(point_count)) {
        omitted = metric_found_path_count - goal_output.path_results.size();
        info.limits_reached =
          static_cast<uint16_t>(info.limits_reached | PlanningResultInfo::LIMIT_MAX_PATH_POINTS);
        break;
      }
      size_t path_bytes = 0;
      size_t topology_path_bytes = 0;
      size_t candidate_response_bytes = response_bytes;
      const auto bytes_fit = [&](size_t serialized_path_points, size_t& candidate_bytes) {
        candidate_bytes = response_bytes;
        return estimatePathResponseBytes(
                 serialized_path_points, grid.header.frame_id.size(), path_bytes) &&
               estimatePathResponseBytes(
                 topology_point_count, grid.header.frame_id.size(), topology_path_bytes) &&
               checkedAdd(candidate_bytes, path_bytes, configuration.planning.max_response_bytes) &&
               checkedAdd(
                 candidate_bytes, topology_path_bytes, configuration.planning.max_response_bytes);
      };
      bool response_bytes_fit = bytes_fit(point_count, candidate_response_bytes);
      if (!response_bytes_fit && !topology_only) {
        topology_only = true;
        point_count = topology_point_count;
        response_bytes_fit = bytes_fit(point_count, candidate_response_bytes);
      }
      if (!response_bytes_fit) {
        omitted = metric_found_path_count - goal_output.path_results.size();
        info.limits_reached =
          static_cast<uint16_t>(info.limits_reached | PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
        break;
      }
      if (topology_only) {
        path_result.path = path_result.topology_path;
      } else if (!buildPathMsg(solution,
                               work_map,
                               grid.header.frame_id,
                               configuration.planning.max_path_points - emitted_path_points -
                                 topology_point_count,
                               path_result.path,
                               path_error)) {
        ++omitted;
        continue;
      }
      const double core_metric_cost = gridCostToMetric(solution.path_cost_, work_map.resolution);
      if (!finalizeCostBoundedPublishedPathResult(
            path_result, core_metric_cost, request.max_path_lengths[i], stop_requested)) {
        ++omitted;
        continue;
      }
      // As in the single-goal path, admission accounting must follow any
      // dense-to-topology fallback performed by the final certificate.
      if (!bytes_fit(path_result.path.poses.size(), candidate_response_bytes)) {
        ++omitted;
        info.limits_reached =
          static_cast<uint16_t>(info.limits_reached | PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
        continue;
      }
      const size_t serialized_path_points =
        path_result.path.poses.size() + path_result.topology_path.poses.size();
      response_bytes = candidate_response_bytes;
      goal_output.path_results.emplace_back(std::move(path_result));
      goal_emitted_solution_indices.emplace_back(solution_index);
      emitted_path_points += serialized_path_points;
      ++total_returned_paths;
    }
    if (!stableSortPublishedPathsWithSources(goal_output.path_results,
                                             goal_emitted_solution_indices)) {
      throw std::logic_error("published goal path/source index cardinality mismatch");
    }
    for (const size_t solution_index : goal_emitted_solution_indices)
      emitted_solutions.emplace_back(i, solution_index);
    info.returned_path_count = boundedUint32(goal_output.path_results.size());
    info.output_complete =
      omitted == 0 && goal_output.path_results.size() == metric_found_path_count;
    if (!info.output_complete) {
      markPathOutputIncomplete(info);
      all_output_complete = false;
    }
    info.request_satisfied = info.cost_bound_exhausted && info.output_complete;
    aggregate.limits_reached =
      static_cast<uint16_t>(aggregate.limits_reached | info.limits_reached);
    goal_output.success = !goal_output.path_results.empty();
    if (info.search_complete && metric_found_path_count == 0)
      info.status = PlanningResultInfo::STATUS_NO_PATH;
    if (goal_output.success)
      ++goals_with_paths;
    all_requests_satisfied = all_requests_satisfied && info.request_satisfied;
    response.goal_results.emplace_back(std::move(goal_output));
  }

  const bool complete_goal_result_set = result.goal_results.size() == request.goals.size() &&
                                        response.goal_results.size() == request.goals.size();
  if (!complete_goal_result_set) {
    all_output_complete = false;
    all_requests_satisfied = false;
  }
  aggregate.returned_goal_count = boundedUint32(response.goal_results.size());
  aggregate.completed_goal_count = boundedUint32(completed_goals);
  aggregate.goals_with_paths = boundedUint32(goals_with_paths);
  aggregate.found_path_count = boundedUint32(total_found_paths);
  aggregate.returned_path_count = boundedUint32(total_returned_paths);
  aggregate.search_complete = complete_goal_result_set && all_search_complete;
  aggregate.output_complete = complete_goal_result_set && all_output_complete;
  aggregate.request_satisfied = complete_goal_result_set && all_requests_satisfied &&
                                aggregate.search_complete && aggregate.output_complete;
  // A complete bounded enumeration is a successful aggregate operation even
  // when one or every goal has zero admissible paths. Reachability remains a
  // per-goal payload fact rather than changing the Action transport outcome.
  response.success = aggregate.request_satisfied;
  if (aggregate.search_complete && total_found_paths == 0)
    aggregate.status = PlanningResultInfo::STATUS_NO_PATH;
  if (!all_output_complete && (aggregate.status == PlanningResultInfo::STATUS_COMPLETE ||
                               aggregate.status == PlanningResultInfo::STATUS_NO_PATH)) {
    aggregate.status = PlanningResultInfo::STATUS_PARTIAL_OUTPUT;
  }

  const size_t response_budget_remaining =
    configuration.planning.max_response_bytes > response_bytes
      ? configuration.planning.max_response_bytes - response_bytes
      : 0;
  const size_t debug_node_count = request.include_debug
                                    ? std::min({nodes.size(),
                                                configuration.planning.max_debug_nodes,
                                                response_budget_remaining / kDebugNodeBytes})
                                    : 0u;
  response.debug_nodes.reserve(debug_node_count);
  for (size_t i = 0; i < debug_node_count; ++i) {
    if (stop_requested && stop_requested()) {
      mark_goal_set_cancelled();
      return;
    }
    raystar_interfaces::msg::DebugNode debug_node;
    debug_node.index = nodes[i].index();
    debug_node.g_cost = nodes[i].gCost() * resolution;
    debug_node.f_cost = (nodes[i].gCost() + nodes[i].hCost()) * resolution;
    response.debug_nodes.emplace_back(std::move(debug_node));
  }
  response_bytes += debug_node_count * kDebugNodeBytes;
  aggregate.debug_output_complete =
    !request.include_debug || debug_node_count == result.expanded_nodes;
  for (auto& goal_result : response.goal_results)
    goal_result.result_info.debug_output_complete = aggregate.debug_output_complete;

  if (response.success && result.polymap) {
    std::vector<PathSolution> visualization_solutions;
    visualization_solutions.reserve(emitted_solutions.size());
    for (const auto& emitted : emitted_solutions) {
      visualization_solutions.emplace_back(
        result.goal_results[emitted.first].path_solutions[emitted.second]);
    }
    last_frame_id_ = grid.header.frame_id;
    publishNonHomotopicPaths(visualization_solutions,
                             work_map,
                             grid.header.frame_id,
                             configuration.planning.max_path_points,
                             configuration.planning.max_response_bytes);
    if (request.include_debug) {
      publishPolyObstacles(
        *result.polymap, work_map, grid.header.frame_id, configuration.planning.max_response_bytes);
      publishCDT(
        *result.polymap, work_map, grid.header.frame_id, configuration.planning.max_response_bytes);
      publishDebugTree(nodes,
                       work_map,
                       grid.header.frame_id,
                       configuration.planning.max_debug_nodes,
                       configuration.planning.max_response_bytes);
    }
  }
} catch (const std::exception& exception) {
  std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
  clearVisualizationsLocked();
  response.success = false;
  response.result_info.map_id = map_id;
  response.result_info.status = PlanningResultInfo::STATUS_FAILED;
  response.goal_results.clear();
  response.debug_nodes.clear();
  setBoundedExceptionMessage(response.message, exception.what());
} catch (...) {
  std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
  clearVisualizationsLocked();
  response.success = false;
  response.result_info.map_id = map_id;
  response.result_info.status = PlanningResultInfo::STATUS_FAILED;
  response.goal_results.clear();
  response.debug_nodes.clear();
  response.message = "Raystar goal-set request failed with an unknown exception";
}

void RaystarNode::executeTransitionPlanning(
  const TransitionAction::Goal& request,
  TransitionAction::Result& response,
  const nav_msgs::msg::OccupancyGrid& grid,
  const raystar_interfaces::MapId& map_id,
  const StopPredicate& stop_requested,
  const TransitionProgressCallback& progress_callback) noexcept try {
  response = TransitionAction::Result{};
  response.map_id = map_id;
  response.requested_transition_count = boundedUint32(request.transition_pairs.size());
  TransitionProgressReporter progress(request.transition_pairs.size(), progress_callback);
  progress.publishStage("validating transition request");

  RequestConfiguration configuration;
  std::string error;
  if (!loadRequestConfiguration(configuration, error)) {
    response.status = TransitionAction::Result::STATUS_FAILED;
    response.message = std::move(error);
    return;
  }
  response.environment_id = raystar_interfaces::computeEnvironmentId(
    grid, configuration.occupied_threshold, request.allow_unknown);
  if (!raystar_interfaces::isZeroEnvironmentId(request.expected_environment_id) &&
      !raystar_interfaces::environmentIdsEqual(request.expected_environment_id,
                                               response.environment_id)) {
    response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
    response.message =
      "expected_environment_id does not match the cached map, occupancy policy, "
      "or Raystar planning semantics";
    return;
  }
  const bool references_must_be_taut =
    request.reference_path_policy == TransitionAction::Goal::REFERENCE_PATHS_MUST_BE_TAUT;
  if (!references_must_be_taut &&
      request.reference_path_policy != TransitionAction::Goal::REFERENCE_PATHS_MAY_BE_UNTAUT) {
    response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
    response.message =
      "Invalid reference_path_policy: expected REFERENCE_PATHS_MUST_BE_TAUT or "
      "REFERENCE_PATHS_MAY_BE_UNTAUT";
    return;
  }
  const auto planning_deadline =
    std::chrono::steady_clock::now() + configuration.planning.planning_timeout;
  bool planning_timed_out = false;
  const StopPredicate effective_stop_requested = [&]() {
    if (stop_requested && stop_requested())
      return true;
    if (std::chrono::steady_clock::now() >= planning_deadline) {
      planning_timed_out = true;
      return true;
    }
    return false;
  };
  const auto stopped_status = [&]() {
    return planning_timed_out ? TransitionAction::Result::STATUS_TIMEOUT
                              : TransitionAction::Result::STATUS_CANCELLED;
  };
  const auto stopped_reason = [&]() {
    return planning_timed_out ? std::string("timed out") : std::string("was cancelled");
  };
  if (effective_stop_requested()) {
    response.status = stopped_status();
    response.message = "UPS transition construction " + stopped_reason() + " before validation";
    return;
  }
  if (request.rooted_reference_paths.empty()) {
    response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
    response.message = "At least one rooted reference is required";
    return;
  }
  if (request.rooted_reference_paths.size() > configuration.planning.max_transition_references) {
    response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
    response.message = "Rooted reference count exceeds max_transition_references=" +
                       std::to_string(configuration.planning.max_transition_references);
    return;
  }
  if (request.transition_pairs.empty()) {
    response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
    response.message = "At least one directed transition pair is required";
    return;
  }
  if (request.transition_pairs.size() > configuration.planning.max_transition_pairs) {
    response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
    response.message = "Transition pair count exceeds max_transition_pairs=" +
                       std::to_string(configuration.planning.max_transition_pairs);
    return;
  }
  for (size_t index = 0; index < request.transition_pairs.size(); ++index) {
    const auto& pair = request.transition_pairs[index];
    if (pair.from_reference >= request.rooted_reference_paths.size() ||
        pair.to_reference >= request.rooted_reference_paths.size()) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Transition pair " + std::to_string(index) +
                         " contains an out-of-range configuration index";
      return;
    }
  }

  progress.publishStage("preparing transition environment");
  GridMap work_map;
  if (!occupancyGridToBinaryMap(
        grid, request.allow_unknown, configuration, effective_stop_requested, work_map, error)) {
    response.status = effective_stop_requested() ? stopped_status()
                                                 : TransitionAction::Result::STATUS_INVALID_REQUEST;
    response.message = std::move(error);
    return;
  }
  // Direct Polymap construction, unlike RaystarCore::plan(), does not add
  // the planner's mandatory occupied outer ring. Apply the same geometric
  // boundary contract before extracting contours and building the CDT.
  for (unsigned int x = 0; x < work_map.width; ++x) {
    work_map.data[x] = 1;
    work_map.data[(static_cast<size_t>(work_map.height) - 1) * work_map.width + x] = 1;
  }
  for (unsigned int y = 0; y < work_map.height; ++y) {
    work_map.data[static_cast<size_t>(y) * work_map.width] = 1;
    work_map.data[static_cast<size_t>(y) * work_map.width + work_map.width - 1] = 1;
  }
  if (grid.header.frame_id.size() > kMaxFrameIdBytes) {
    response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
    response.message = "Invalid map: frame_id must be at most 256 bytes";
    return;
  }

  const auto waypoint_to_grid = [&work_map](double wx, double wy, Point2d& point) {
    if (!hasValidWorldTransform(work_map) || !std::isfinite(wx) || !std::isfinite(wy))
      return false;
    double gx = (wx - work_map.origin_x) / static_cast<double>(work_map.resolution);
    double gy = (wy - work_map.origin_y) / static_cast<double>(work_map.resolution);
    gx = canonicalizeWorldGridCoordinate(
      wx, work_map.origin_x, static_cast<double>(work_map.resolution), work_map.width, gx);
    gy = canonicalizeWorldGridCoordinate(
      wy, work_map.origin_y, static_cast<double>(work_map.resolution), work_map.height, gy);
    if (!std::isfinite(gx) || !std::isfinite(gy) || gx < 0.0 || gy < 0.0 ||
        gx > static_cast<double>(work_map.width) || gy > static_cast<double>(work_map.height)) {
      return false;
    }
    point = {gx, gy};
    return true;
  };

  std::vector<std::vector<Point2d>> configurations;
  configurations.reserve(request.rooted_reference_paths.size());
  std::vector<PolymapEndpoint> endpoints;
  endpoints.reserve(request.rooted_reference_paths.size());
  size_t input_point_count = 0;
  std::optional<ContinuousGridPoint> root;
  for (size_t configuration_index = 0; configuration_index < request.rooted_reference_paths.size();
       ++configuration_index) {
    const auto& input_path = request.rooted_reference_paths[configuration_index];
    if (input_path.header.frame_id != grid.header.frame_id ||
        input_path.header.frame_id.size() > kMaxFrameIdBytes) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Rooted reference " + std::to_string(configuration_index) +
                         " Path frame_id does not match the cached map";
      return;
    }
    if (input_path.poses.empty()) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message =
        "Rooted reference " + std::to_string(configuration_index) + " is empty";
      return;
    }
    if (!canAppendCount(
          input_point_count, input_path.poses.size(), configuration.planning.max_path_points)) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Input rooted reference paths exceed max_path_points=" +
                         std::to_string(configuration.planning.max_path_points);
      return;
    }
    input_point_count += input_path.poses.size();
    std::vector<Point2d> path;
    path.reserve(input_path.poses.size());
    for (size_t point_index = 0; point_index < input_path.poses.size(); ++point_index) {
      const auto& pose = input_path.poses[point_index];
      std::string pose_error;
      if (!validatePlanarPose(
            pose, "Rooted reference waypoint", grid.header.frame_id, pose_error)) {
        response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
        response.message = "Configuration " + std::to_string(configuration_index) + " waypoint " +
                           std::to_string(point_index) + ": " + pose_error;
        return;
      }
      Point2d point;
      if (!waypoint_to_grid(pose.pose.position.x, pose.pose.position.y, point)) {
        response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
        response.message = "Configuration " + std::to_string(configuration_index) + " waypoint " +
                           std::to_string(point_index) + " lies outside the map geometry";
        return;
      }
      if (path.empty() || path.back() != point)
        path.push_back(point);
    }

    ContinuousGridPoint reference_root;
    const auto& first_pose = input_path.poses.front().pose.position;
    if (!worldToContinuousMap(work_map, first_pose.x, first_pose.y, reference_root)) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Rooted reference " + std::to_string(configuration_index) +
                         " root is not a strict map-interior point";
      return;
    }
    if (!root) {
      root = reference_root;
    } else if (root->x != reference_root.x || root->y != reference_root.y) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "All rooted references must have one identical root point";
      return;
    }
    ContinuousGridPoint endpoint;
    const auto& last_pose = input_path.poses.back().pose.position;
    if (!worldToContinuousMap(work_map, last_pose.x, last_pose.y, endpoint)) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Rooted reference " + std::to_string(configuration_index) +
                         " endpoint is not a strict map-interior point";
      return;
    }
    const PolymapEndpoint polymap_endpoint{static_cast<int>(endpoint.cell_x),
                                           static_cast<int>(endpoint.cell_y),
                                           {endpoint.x, endpoint.y}};
    const auto duplicate_endpoint =
      std::find_if(endpoints.begin(), endpoints.end(), [&polymap_endpoint](const auto& existing) {
        return existing.cell_x == polymap_endpoint.cell_x &&
               existing.cell_y == polymap_endpoint.cell_y &&
               existing.position == polymap_endpoint.position;
      });
    if (duplicate_endpoint == endpoints.end())
      endpoints.push_back(polymap_endpoint);
    configurations.emplace_back(std::move(path));
  }

  const StopToken stop_token(effective_stop_requested);
  const PolymapEndpoint root_endpoint{
    static_cast<int>(root->cell_x), static_cast<int>(root->cell_y), Point2d{root->x, root->y}};
  // Only transition construction writes this cache. Planner Polymaps use
  // simplified contours and must never satisfy a raw-contour UPS lookup.
  std::shared_ptr<const Polymap> polymap_owner = findCachedTransitionEnvironment(
    grid, map_id, request.allow_unknown, configuration, root_endpoint, endpoints);
  if (polymap_owner) {
    RCLCPP_DEBUG(get_logger(),
                 "Reusing the completed free-triangle environment for %zu UPS transition(s)",
                 request.transition_pairs.size());
  } else {
    auto polymap_result = Polymap::createForReferenceShortening(work_map,
                                                                 root_endpoint.cell_x,
                                                                 root_endpoint.cell_y,
                                                                 root_endpoint.position,
                                                                 endpoints,
                                                                 stop_token,
                                                                 configuration.planning);
    if (polymap_result.status == PolymapCreateStatus::stopped) {
      response.status = stopped_status();
      response.message =
        "UPS transition construction " + stopped_reason() + " while building the map";
      return;
    }
    if (!polymap_result) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = polymap_result.error.empty()
                           ? "Could not build the free-triangle environment"
                           : polymap_result.error.substr(0, kMaxDiagnosticBytes);
      return;
    }
    polymap_owner = std::make_shared<Polymap>(std::move(*polymap_result.value));
    cacheCompletedTransitionEnvironment(
      grid, map_id, request.allow_unknown, configuration, root_endpoint, endpoints, polymap_owner);
  }
  const Polymap& polymap = *polymap_owner;

  // Validate every complete configuration before any pairwise shortening.
  // shortenWithinHomotopy() deliberately removes an exact common prefix from
  // alpha_a and alpha_b. Without this admission pass, two caller-supplied
  // references could share an obstacle-crossing prefix that never appears in
  // the composed alpha_a^{-1} * alpha_b path and would therefore evade the
  // pairwise corridor trace.
  progress.publishStage("validating rooted references");
  for (size_t configuration_index = 0; configuration_index < configurations.size();
       ++configuration_index) {
    if (stop_token.poll()) {
      response.status = stopped_status();
      response.message = "UPS transition construction " + stopped_reason() +
                         " while validating rooted reference " +
                         std::to_string(configuration_index);
      return;
    }
    const auto validation =
      polymap.shortenPathWithinHomotopy(configurations[configuration_index], stop_token);
    if (validation.status == HomotopyShorteningStatus::stopped) {
      response.status = stopped_status();
      response.message = "UPS transition construction " + stopped_reason() +
                         " while validating rooted reference " +
                         std::to_string(configuration_index);
      return;
    }
    if (!validation) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Rooted reference " + std::to_string(configuration_index) +
                         " is not a collision-free reference in the cached map";
      if (!validation.message.empty())
        response.message += ": " + validation.message.substr(0, kMaxDiagnosticBytes);
      response.message.resize(std::min(response.message.size(), kMaxDiagnosticBytes));
      return;
    }
    if (references_must_be_taut) {
      double input_cost = 0.0;
      for (size_t point_index = 1; point_index < configurations[configuration_index].size();
           ++point_index) {
        const auto& previous = configurations[configuration_index][point_index - 1];
        const auto& current = configurations[configuration_index][point_index];
        input_cost += std::hypot(current.first - previous.first, current.second - previous.second);
      }
      const double taut_tolerance = 1.0e-10 * std::max({1.0, input_cost, validation.path_cost});
      if (!std::isfinite(input_cost) ||
          std::abs(input_cost - validation.path_cost) > taut_tolerance) {
        response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
        response.message = "Rooted reference " + std::to_string(configuration_index) +
                           " is not a locally shortest (taut) reference";
        return;
      }
    }
  }

  response.transitions.reserve(request.transition_pairs.size());
  size_t output_point_count = 0;
  size_t response_bytes = kResponseBaseBytes;
  bool all_successful = true;
  size_t failed_transition_count = 0;
  size_t first_failed_index = 0;
  uint8_t first_failed_status = 0;
  std::string first_failed_message;
  progress.publishStage("shortening transition pairs");
  for (size_t index = 0; index < request.transition_pairs.size(); ++index) {
    if (stop_token.poll()) {
      response.status = stopped_status();
      response.message = "UPS transition construction " + stopped_reason() + " after " +
                         std::to_string(response.transitions.size()) + " transitions";
      return;
    }
    const auto& requested_pair = request.transition_pairs[index];
    const auto shortening =
      RaystarCore::shortenWithinHomotopy(polymap,
                                         configurations[requested_pair.from_reference],
                                         configurations[requested_pair.to_reference],
                                         stop_token);
    if (shortening.status == HomotopyShorteningStatus::stopped) {
      response.status = stopped_status();
      response.message = "UPS transition construction " + stopped_reason() +
                         " while shortening pair " + std::to_string(index);
      return;
    }

    raystar_interfaces::msg::HomotopyTransitionResult output;
    output.pair = requested_pair;
    switch (shortening.status) {
    case HomotopyShorteningStatus::success:
      output.status = output.STATUS_SUCCESS;
      break;
    case HomotopyShorteningStatus::invalid_reference:
      output.status = output.STATUS_INVALID_REFERENCE;
      break;
    case HomotopyShorteningStatus::no_corridor:
      output.status = output.STATUS_NO_CORRIDOR;
      break;
    case HomotopyShorteningStatus::stopped:
      output.status = output.STATUS_STOPPED;
      break;
    case HomotopyShorteningStatus::failure:
      output.status = output.STATUS_FAILURE;
      break;
    }
    output.collision_free = shortening.collision_free;
    output.homotopy_preserved = shortening.homotopy_preserved;
    output.locally_shortest = shortening.locally_shortest;
    output.triangle_occurrences = shortening.corridor.triangle_occurrences;
    output.message = shortening.message.substr(0, kMaxDiagnosticBytes);
    if (!shortening) {
      if (failed_transition_count == 0) {
        first_failed_index = index;
        first_failed_status = output.status;
        first_failed_message = output.message;
      }
      ++failed_transition_count;
    }
    output.path.header.stamp = now();
    output.path.header.frame_id = grid.header.frame_id;

    if (!canAppendCount(
          output_point_count, shortening.path.size(), configuration.planning.max_path_points)) {
      response.status = TransitionAction::Result::STATUS_FAILED;
      response.message = "UPS output exceeds max_path_points=" +
                         std::to_string(configuration.planning.max_path_points);
      return;
    }
    size_t path_bytes = 0;
    const bool triangle_bytes_valid =
      output.triangle_occurrences.size() <= std::numeric_limits<size_t>::max() / sizeof(uint32_t);
    if (!estimatePathResponseBytes(
          shortening.path.size(), grid.header.frame_id.size(), path_bytes) ||
        !triangle_bytes_valid ||
        !checkedAdd(response_bytes, path_bytes, configuration.planning.max_response_bytes) ||
        !checkedAdd(
          response_bytes, output.message.size(), configuration.planning.max_response_bytes) ||
        !checkedAdd(response_bytes,
                    output.triangle_occurrences.size() * sizeof(uint32_t),
                    configuration.planning.max_response_bytes)) {
      response.status = TransitionAction::Result::STATUS_FAILED;
      response.message = "UPS output exceeds max_response_bytes=" +
                         std::to_string(configuration.planning.max_response_bytes);
      return;
    }
    output.path.poses.reserve(shortening.path.size());
    for (const auto& point : shortening.path) {
      geometry_msgs::msg::PoseStamped pose;
      pose.header = output.path.header;
      pose.pose.orientation.w = 1.0;
      if (!continuousGridToWorld(work_map, point, pose.pose.position.x, pose.pose.position.y)) {
        response.status = TransitionAction::Result::STATUS_FAILED;
        response.message = "Could not convert a UPS output waypoint to world coordinates";
        return;
      }
      output.path.poses.emplace_back(std::move(pose));
    }
    if (shortening) {
      double rounded_world_length = 0.0;
      double certified_world_length = 0.0;
      const auto certificate_stop = [&stop_token]() { return stop_token.poll(); };
      if (output.path.poses.size() > 1 &&
          !publishedPathLength(
            output.path, rounded_world_length, certified_world_length, certificate_stop)) {
        const bool stopped = stop_token.poll();
        response.status = stopped ? stopped_status() : TransitionAction::Result::STATUS_FAILED;
        response.message = stopped ? "UPS transition length certification " + stopped_reason()
                                   : "Could not certify the serialized UPS path length";
        return;
      }
      const double nominal_metric_length =
        gridCostToMetric(shortening.path_cost, work_map.resolution);
      if (!std::isfinite(nominal_metric_length) ||
          !publishedLengthsMatch(nominal_metric_length, certified_world_length) ||
          (output.path.poses.size() == 1 && nominal_metric_length != 0.0)) {
        response.status = TransitionAction::Result::STATUS_FAILED;
        response.message =
          "Serialized UPS path length is inconsistent with the certified Core path";
        return;
      }
      output.path_length = certified_world_length;
    }
    output_point_count += shortening.path.size();
    all_successful = all_successful && static_cast<bool>(shortening);
    response.transitions.emplace_back(std::move(output));
    response.completed_transition_count = boundedUint32(response.transitions.size());
    progress.publishPairCompletion(response.transitions.size(), "shortening transition pairs");
  }

  response.success = all_successful;
  response.status = all_successful ? TransitionAction::Result::STATUS_COMPLETE
                                   : TransitionAction::Result::STATUS_FAILED;
  if (all_successful) {
    response.message =
      "All requested UPS transitions are collision-free, homotopy-preserving, "
      "and locally shortest";
  } else {
    std::ostringstream diagnostic;
    diagnostic << failed_transition_count << " of " << request.transition_pairs.size()
               << " UPS transitions could not be certified; first failure index="
               << first_failed_index << " status=" << static_cast<unsigned int>(first_failed_status)
               << ": " << first_failed_message;
    response.message = diagnostic.str().substr(0, kMaxDiagnosticBytes);
  }
  progress.publishFinal("transition batch complete");
} catch (const std::exception& exception) {
  response.success = false;
  response.map_id = map_id;
  response.status = TransitionAction::Result::STATUS_FAILED;
  response.transitions.clear();
  response.completed_transition_count = 0;
  setBoundedExceptionMessage(response.message, exception.what());
} catch (...) {
  response.success = false;
  response.map_id = map_id;
  response.status = TransitionAction::Result::STATUS_FAILED;
  response.transitions.clear();
  response.completed_transition_count = 0;
  response.message = "Raystar UPS request failed with an unknown exception";
}


// executePlanning is defined in this translation unit and called from
// raystar_node.cpp (GetRaystarPaths service) and raystar_node_actions.cpp
// (PlanRaystarPaths action).  Instantiate both request/response pairs here so
// the template definition can stay out of raystar_node.h.
template void RaystarNode::executePlanning<raystar_interfaces::srv::GetRaystarPaths::Request,
                                           raystar_interfaces::srv::GetRaystarPaths::Response>(
  const raystar_interfaces::srv::GetRaystarPaths::Request& request,
  raystar_interfaces::srv::GetRaystarPaths::Response& response,
  const nav_msgs::msg::OccupancyGrid& grid,
  const raystar_interfaces::MapId& map_id,
  const StopPredicate& stop_requested) noexcept;

template void RaystarNode::executePlanning<RaystarNode::PlanAction::Goal,
                                           RaystarNode::PlanAction::Result>(
  const RaystarNode::PlanAction::Goal& request,
  RaystarNode::PlanAction::Result& response,
  const nav_msgs::msg::OccupancyGrid& grid,
  const raystar_interfaces::MapId& map_id,
  const StopPredicate& stop_requested) noexcept;

}  // namespace raystar
