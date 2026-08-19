#ifndef RAYSTAR_CORE_DETAIL_H_
#define RAYSTAR_CORE_DETAIL_H_

// Package-private helpers shared by the RaystarCore translation units
// (raystar_core.cpp, raystar_core_planning.cpp).  Formerly an anonymous
// namespace inside the single raystar_core.cpp.

#include <raystar/raystar_core.h>
#include <chrono>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <limits>
#include <exception>
#include <optional>

#include "exact_point_location.h"
#include "conservative_path_length.h"
#include "visibility_validation.h"


namespace raystar {

namespace core_impl {

inline detail::VisibilityBoundaryMode visibilityBoundaryModeForNode(const Polymap& polymap,
                                                                    const Node& node) {
  if (node.isContinuousRoot())
    return detail::VisibilityBoundaryMode::closed_cycle;
  if (node.parentIndex() >= 0)
    return detail::VisibilityBoundaryMode::open_sector;
  if (node.hasValidSeed() &&
      polymap.isValidTopology(polymap.locateVertex(node.seed().first, node.seed().second))) {
    return detail::VisibilityBoundaryMode::open_sector;
  }
  return detail::VisibilityBoundaryMode::closed_cycle;
}


inline bool isRepresentableAsInt(double value) {
  return std::isfinite(value) && value >= static_cast<double>(std::numeric_limits<int>::lowest()) &&
         value <= static_cast<double>(std::numeric_limits<int>::max());
}

inline int safeIntCast(double value) {
  return isRepresentableAsInt(value) ? static_cast<int>(value) : 0;
}

inline int safeRoundedIntCast(double value) {
  if (!isRepresentableAsInt(value))
    return 0;
  const double rounded = std::round(value);
  return isRepresentableAsInt(rounded) ? static_cast<int>(rounded) : 0;
}

inline bool addCertifiedSegment(ConservativeBinary64PathLength& certificate,
                         const Point2d& first,
                         const Point2d& second) {
  return certificate.addSegment(first.first, first.second, second.first, second.second);
}

inline bool buildPolylineLengthCertificate(const Point2d& start,
                                    const std::vector<std::pair<int, int>>& turning_points,
                                    const Point2d& goal,
                                    ConservativeBinary64PathLength& certificate,
                                    const StopToken* stop_token = nullptr) {
  Point2d previous = start;
  for (const auto& turning_point : turning_points) {
    if (stop_token && stop_token->poll())
      return false;
    const Point2d current = {static_cast<double>(turning_point.first),
                             static_cast<double>(turning_point.second)};
    if (!addCertifiedSegment(certificate, previous, current))
      return false;
    previous = current;
  }
  if (stop_token && stop_token->poll())
    return false;
  return addCertifiedSegment(certificate, previous, goal);
}

inline bool buildCandidateLengthCertificate(const Point2d& start,
                                     const std::vector<Node>& nodes,
                                     const Candidate& candidate,
                                     const Point2d& goal,
                                     ConservativeBinary64PathLength& certificate,
                                     const StopToken* stop_token = nullptr) {
  if (stop_token && stop_token->poll())
    return false;
  if (candidate.nodeIndex() < 0)
    return addCertifiedSegment(certificate, start, goal);
  const size_t node_index = static_cast<size_t>(candidate.nodeIndex());
  if (node_index >= nodes.size())
    return false;
  const auto& node = nodes[node_index];
  const int child_index_value = candidate.childIndex();
  if (child_index_value < 0 || static_cast<size_t>(child_index_value) >= node.children().size()) {
    return false;
  }

  Point2d previous = start;
  for (const auto& turning_point : node.localShortestPath()) {
    if (stop_token && stop_token->poll())
      return false;
    const Point2d current = {static_cast<double>(turning_point.first),
                             static_cast<double>(turning_point.second)};
    if (!addCertifiedSegment(certificate, previous, current))
      return false;
    previous = current;
  }
  if (stop_token && stop_token->poll())
    return false;
  const Point2d& child =
    node.children()[static_cast<size_t>(child_index_value)].endpoint().position;
  return addCertifiedSegment(certificate, previous, child) &&
         addCertifiedSegment(certificate, child, goal);
}

// Validate the map shape before any caller-owned GridMap is copied into the
// planner's working map.  Besides making the multiplication overflow-safe,
// this is the single admission point for the configurable cell budget.  The
// signed-indexing check remains an independent invariant because several
// legacy geometry tables use int coordinates even when the configured budget
// is larger.
inline bool validateMapShapeAndAdmission(const GridMap& grid_map,
                                  const PlanningLimits& limits,
                                  size_t& width,
                                  size_t& height,
                                  size_t& cell_count,
                                  std::string& error) {
  width = static_cast<size_t>(grid_map.width);
  height = static_cast<size_t>(grid_map.height);
  MapResourceEstimate estimate;
  if (!validateMapResourceBudget(width, height, grid_map.data.size(), limits, estimate, error)) {
    return false;
  }
  cell_count = estimate.cell_count;
  return true;
}

inline bool validateNonMapPlanningLimits(const PlanningLimits& limits, std::string& error) {
  const size_t max_int = static_cast<size_t>(std::numeric_limits<int>::max());
  if (limits.max_k <= 0) {
    error = "Invalid planning limits: max_k must be greater than zero";
    return false;
  }
  if (limits.max_nodes == 0 || limits.max_nodes > max_int) {
    error = "Invalid planning limits: max_nodes must be between 1 and " + std::to_string(max_int);
    return false;
  }
  if (limits.max_cost_bounded_paths == 0 || limits.max_cost_bounded_paths > max_int) {
    error = "Invalid planning limits: max_cost_bounded_paths must be between 1 and " +
            std::to_string(max_int);
    return false;
  }
  if (limits.max_multi_goal_count == 0 || limits.max_multi_goal_count > max_int) {
    error = "Invalid planning limits: max_multi_goal_count must be between 1 and " +
            std::to_string(max_int);
    return false;
  }
  if (limits.planning_timeout.count() < 0) {
    error = "Invalid planning limits: planning timeout must not be negative";
    return false;
  }
  if (limits.max_path_points < 2) {
    error = "Invalid planning limits: max_path_points must be at least 2";
    return false;
  }
  // Zero is an intentional opt-out for the ROS debug arrays/markers.  It
  // never affects the search-node admission limit itself.
  if (limits.max_debug_nodes > max_int) {
    error =
      "Invalid planning limits: max_debug_nodes must be between 0 and " + std::to_string(max_int);
    return false;
  }
  if (limits.max_response_bytes == 0) {
    error = "Invalid planning limits: max_response_bytes must be greater than zero";
    return false;
  }
  return true;
}

inline bool validateSearchObjective(const SearchObjective& objective,
                             const PlanningLimits& limits,
                             std::string& error) {
  switch (objective.mode) {
  case SearchMode::top_k:
    if (objective.k <= 0) {
      error = "Invalid K: must be greater than zero";
      return false;
    }
    if (objective.k > limits.max_k) {
      error = "Invalid K: requested " + std::to_string(objective.k) +
              " exceeds max_k=" + std::to_string(limits.max_k);
      return false;
    }
    if (objective.max_path_cost != 0.0) {
      error = "Invalid top-K objective: max_path_cost must be zero";
      return false;
    }
    if (objective.path_admission) {
      error = "Invalid top-K objective: bounded path admission must be empty";
      return false;
    }
    break;
  case SearchMode::all_within_cost:
    if (objective.k != 0) {
      error = "Invalid cost-bounded objective: K must be zero";
      return false;
    }
    if (!std::isfinite(objective.max_path_cost) || objective.max_path_cost < 0.0) {
      error = "Invalid cost-bounded objective: max_path_cost must be finite and non-negative";
      return false;
    }
    break;
  default:
    error = "Invalid search objective mode";
    return false;
  }
  error.clear();
  return true;
}

inline std::string boundedPlanningExceptionMessage(const char* detail) {
  constexpr size_t kMaxMessageBytes = 1024;
  constexpr const char* prefix = "Planning failed with exception: ";
  std::string message(prefix);
  if (message.size() >= kMaxMessageBytes || detail == nullptr) {
    message.resize(std::min(message.size(), kMaxMessageBytes));
    return message;
  }
  const size_t remaining = kMaxMessageBytes - message.size();
  size_t detail_bytes = 0;
  while (detail_bytes < remaining && detail[detail_bytes] != '\0') ++detail_bytes;
  message.append(detail, detail_bytes);
  return message;
}

inline bool validatePlanRequest(const GridMap& grid_map,
                         int start_x,
                         int start_y,
                         int goal_x,
                         int goal_y,
                         const SearchObjective& objective,
                         const PlanningLimits& limits,
                         std::string& error) {
  size_t width = 0;
  size_t height = 0;
  size_t cell_count = 0;
  if (!validateMapShapeAndAdmission(grid_map, limits, width, height, cell_count, error)) {
    return false;
  }
  (void)cell_count;
  if (!validateNonMapPlanningLimits(limits, error))
    return false;

  if (!validateSearchObjective(objective, limits, error))
    return false;

  const auto is_in_bounds = [width, height](int x, int y) {
    return x >= 0 && y >= 0 && static_cast<size_t>(x) < width && static_cast<size_t>(y) < height;
  };

  if (!is_in_bounds(start_x, start_y)) {
    error = "Invalid start: coordinates (" + std::to_string(start_x) + ", " +
            std::to_string(start_y) + ") are outside map bounds";
    return false;
  }
  if (!is_in_bounds(goal_x, goal_y)) {
    error = "Invalid goal: coordinates (" + std::to_string(goal_x) + ", " + std::to_string(goal_y) +
            ") are outside map bounds";
    return false;
  }

  error.clear();
  return true;
}

inline std::string planningLimitMessage(PlanningLimitReached limit,
                                 const PlanningLimits& limits,
                                 size_t path_count) {
  std::string message;
  if (limit == PlanningLimitReached::max_nodes) {
    message = "Planning stopped after reaching max_nodes=" + std::to_string(limits.max_nodes);
  } else if (limit == PlanningLimitReached::max_path_points) {
    message =
      "Planning stopped before exceeding max_path_points=" + std::to_string(limits.max_path_points);
  } else if (limit == PlanningLimitReached::max_paths) {
    message = "Planning stopped after reaching max_cost_bounded_paths=" +
              std::to_string(limits.max_cost_bounded_paths);
  } else if (limit == PlanningLimitReached::timeout) {
    message = "Planning stopped after reaching planning timeout of " +
              std::to_string(limits.planning_timeout.count()) + " ms";
  } else if (limit == PlanningLimitReached::cancelled) {
    message = "Planning was cancelled";
  } else {
    message = "Planning stopped";
  }

  if (path_count == 0) {
    message += limit == PlanningLimitReached::cancelled ? "; no path was found before cancellation"
                                                        : "; no path was found before the limit";
  } else {
    message +=
      "; returning " + std::to_string(path_count) + " path(s), but the search is incomplete";
  }
  return message;
}

inline std::pair<int, int> legacyTopology(const BoundaryEndpoint& endpoint) {
  if (const auto* vertex = std::get_if<ObstacleVertexId>(&endpoint.support))
    return {vertex->obstacle, vertex->vertex};
  return {-1, -1};
}

inline bool projectVisibilityRegion(const VisibilityRegion& visibility_region,
                             std::vector<std::pair<double, double>>& positions,
                             std::vector<std::pair<int, int>>& topology,
                             const StopToken& stop_token) {
  positions.clear();
  topology.clear();
  if (stop_token.poll())
    return false;
  positions.reserve(visibility_region.size());
  topology.reserve(visibility_region.size());
  for (const auto& endpoint : visibility_region) {
    if (stop_token.poll()) {
      positions.clear();
      topology.clear();
      return false;
    }
    positions.emplace_back(endpoint.position);
    topology.emplace_back(legacyTopology(endpoint));
  }
  if (stop_token.poll()) {
    positions.clear();
    topology.clear();
    return false;
  }
  return true;
}

inline BoundaryEndpoint exactEndpoint(const std::pair<int, int>& position,
                               const std::pair<int, int>& topology) {
  BoundaryEndpoint endpoint;
  endpoint.position = {static_cast<double>(position.first), static_cast<double>(position.second)};
  endpoint.exact_position = exact_geometry::Point(position.first, position.second);
  endpoint.support = ObstacleVertexId{topology.first, topology.second};
  return endpoint;
}

inline bool canonicalizeExactObstacleVertices(const Polymap* pMap,
                                       VisibilityRegion& visibility_region,
                                       const StopToken& stop_token) {
  if (!pMap)
    return true;
  if (stop_token.poll())
    return false;
  constexpr double coordinate_tolerance = 1e-9;
  for (auto& endpoint : visibility_region) {
    if (stop_token.poll())
      return false;
    const auto* vertex_id = std::get_if<ObstacleVertexId>(&endpoint.support);
    if (!vertex_id)
      continue;
    const std::pair<int, int> topology = {vertex_id->obstacle, vertex_id->vertex};
    if (!pMap->isValidTopology(topology))
      continue;
    const auto& vertex = pMap->obstacles()[static_cast<size_t>(vertex_id->obstacle)]
                           .ordered_vertices_[static_cast<size_t>(vertex_id->vertex)];
    if (!std::isfinite(endpoint.position.first) || !std::isfinite(endpoint.position.second) ||
        std::abs(endpoint.position.first - static_cast<double>(vertex.first)) >
          coordinate_tolerance ||
        std::abs(endpoint.position.second - static_cast<double>(vertex.second)) >
          coordinate_tolerance) {
      continue;
    }
    endpoint.position = {static_cast<double>(vertex.first), static_cast<double>(vertex.second)};
    endpoint.exact_position = exact_geometry::Point(vertex.first, vertex.second);
  }
  return true;
}

inline VisibilityRegion adaptLegacyVisibilityRegion(
  const Polymap* pMap,
  const std::vector<std::pair<double, double>>& positions,
  const std::vector<std::pair<int, int>>& topology,
  bool& valid,
  std::string& error) {
  VisibilityRegion result;
  valid = true;
  error.clear();

  if (!pMap) {
    valid = false;
    error = "Cannot adapt a legacy visibility region without a Polymap";
    return result;
  }
  if (positions.size() != topology.size()) {
    valid = false;
    error = "Visibility region and topology sizes do not match";
    return result;
  }

  result.reserve(positions.size());
  for (size_t i = 0; i < positions.size(); ++i) {
    BoundaryEndpoint endpoint;
    endpoint.position = positions[i];
    if (!isRepresentableAsInt(endpoint.position.first) ||
        !isRepresentableAsInt(endpoint.position.second)) {
      valid = false;
      error =
        "Visibility region contains a vertex outside the representable int range "
        "at position " +
        std::to_string(i);
      return {};
    }
    endpoint.exact_position =
      exact_geometry::Point(endpoint.position.first, endpoint.position.second);

    const auto& legacy = topology[i];
    if (Polymap::isNoTopology(legacy)) {
      endpoint.support = std::monostate{};
    } else if (!pMap->isValidTopology(legacy)) {
      endpoint.support = ObstacleVertexId{legacy.first, legacy.second};
    } else {
      const auto& obstacle = pMap->obstacles()[static_cast<size_t>(legacy.first)];
      const auto& vertex = obstacle.ordered_vertices_[static_cast<size_t>(legacy.second)];
      constexpr double coordinate_tolerance = 1e-9;
      if (std::abs(endpoint.position.first - static_cast<double>(vertex.first)) <=
            coordinate_tolerance &&
          std::abs(endpoint.position.second - static_cast<double>(vertex.second)) <=
            coordinate_tolerance) {
        endpoint.position = {static_cast<double>(vertex.first), static_cast<double>(vertex.second)};
        endpoint.exact_position = exact_geometry::Point(vertex.first, vertex.second);
        endpoint.support = ObstacleVertexId{legacy.first, legacy.second};
      } else {
        const int next_vertex =
          (legacy.second + 1) % static_cast<int>(obstacle.ordered_vertices_.size());
        endpoint.support = DirectedObstacleEdge{ObstacleVertexId{legacy.first, legacy.second},
                                                ObstacleVertexId{legacy.first, next_vertex}};
      }
    }
    result.emplace_back(std::move(endpoint));
  }
  return result;
}

inline bool validatePlanEndpointShape(const GridMap& grid_map,
                               const PlanEndpoint& endpoint,
                               const char* label,
                               std::string& error) {
  // Check the geometric value before deriving diagnostics from `cell_`.
  // The Point2d overload intentionally uses {-1,-1} as a sentinel when it
  // receives NaN/Inf; reporting an out-of-bounds cell in that case hides the
  // actual input error.
  if (!std::isfinite(endpoint.position_.first) || !std::isfinite(endpoint.position_.second)) {
    error = std::string("Invalid ") + label + ": coordinates must be finite";
    return false;
  }

  const auto& cell = endpoint.cell_;
  if (cell.first < 0 || cell.second < 0 || static_cast<size_t>(cell.first) >= grid_map.width ||
      static_cast<size_t>(cell.second) >= grid_map.height) {
    error = std::string("Invalid ") + label + ": cell coordinates (" + std::to_string(cell.first) +
            ", " + std::to_string(cell.second) + ") are outside map bounds";
    return false;
  }
  const double floor_x = std::floor(endpoint.position_.first);
  const double floor_y = std::floor(endpoint.position_.second);
  if (floor_x != static_cast<double>(cell.first) || floor_y != static_cast<double>(cell.second)) {
    error = std::string("Invalid ") + label + ": continuous position and grid cell disagree";
    return false;
  }
  const auto index =
    static_cast<size_t>(cell.second) * grid_map.width + static_cast<size_t>(cell.first);
  if (index >= grid_map.data.size() || grid_map.data[index] != 0) {
    error = std::string("Invalid ") + label + ": corresponding occupancy cell is occupied";
    return false;
  }
  return true;
}

inline bool validatePlanRequest(const GridMap& grid_map,
                         const PlanEndpoint& start,
                         const PlanEndpoint& goal,
                         const SearchObjective& objective,
                         const PlanningLimits& limits,
                         std::string& error) {
  // Keep all map/limit/K diagnostics byte-for-byte compatible with the legacy
  // request validator, then apply the endpoint-specific contract.
  size_t width = 0;
  size_t height = 0;
  size_t cell_count = 0;
  if (!validateMapShapeAndAdmission(grid_map, limits, width, height, cell_count, error)) {
    return false;
  }
  (void)cell_count;
  if (!validateNonMapPlanningLimits(limits, error))
    return false;
  if (!validateSearchObjective(objective, limits, error))
    return false;
  if (!validatePlanEndpointShape(grid_map, start, "start", error) ||
      !validatePlanEndpointShape(grid_map, goal, "goal", error)) {
    return false;
  }
  error.clear();
  return true;
}

inline bool validateMultiGoalPlanRequest(const GridMap& grid_map,
                                  const PlanEndpoint& start,
                                  const std::vector<CostBoundedGoal>& goals,
                                  const PlanningLimits& limits,
                                  std::string& error) {
  size_t width = 0;
  size_t height = 0;
  size_t cell_count = 0;
  if (!validateMapShapeAndAdmission(grid_map, limits, width, height, cell_count, error))
    return false;
  (void)width;
  (void)height;
  (void)cell_count;
  if (!validateNonMapPlanningLimits(limits, error))
    return false;
  if (goals.empty()) {
    error = "Invalid multi-goal request: at least one goal is required";
    return false;
  }
  if (goals.size() > limits.max_multi_goal_count) {
    error = "Invalid multi-goal request: requested " + std::to_string(goals.size()) +
            " goals exceeds max_multi_goal_count=" + std::to_string(limits.max_multi_goal_count);
    return false;
  }
  if (!validatePlanEndpointShape(grid_map, start, "start", error))
    return false;
  for (size_t i = 0; i < goals.size(); ++i) {
    const auto& goal = goals[i];
    if (!std::isfinite(goal.max_path_cost) || goal.max_path_cost < 0.0) {
      error = "Invalid goal[" + std::to_string(i) +
              "] cost bound: max_path_cost must be finite and non-negative";
      return false;
    }
    const std::string label = "goal[" + std::to_string(i) + "]";
    if (!validatePlanEndpointShape(grid_map, goal.endpoint, label.c_str(), error))
      return false;
  }
  error.clear();
  return true;
}

inline OperationStatus hybridPathSelfCrosses(const Point2d& start,
                                      const std::vector<std::pair<int, int>>& turns,
                                      const Point2d& candidate,
                                      bool& self_crosses,
                                      const StopToken& stop_token) {
  self_crosses = false;
  if (stop_token.poll())
    return OperationStatus::stopped;
  std::vector<exact_geometry::Point> points;
  points.reserve(turns.size() + 2);
  points.emplace_back(start.first, start.second);
  for (const auto& turn : turns) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    points.emplace_back(turn.first, turn.second);
  }
  points.emplace_back(candidate.first, candidate.second);

  const auto& new_point = points.back();
  for (size_t i = 0; i + 1 < points.size() - 1; ++i) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    if (points[i] == new_point) {
      self_crosses = true;
      return OperationStatus::success;
    }
  }
  if (points.size() >= 4) {
    const auto new_segment =
      exact_geometry::Kernel::Segment_2(points[points.size() - 2], points.back());
    for (size_t i = 0; i + 2 < points.size() - 1; ++i) {
      if (stop_token.poll())
        return OperationStatus::stopped;
      const auto old_segment = exact_geometry::Kernel::Segment_2(points[i], points[i + 1]);
      if (CGAL::do_intersect(old_segment, new_segment)) {
        self_crosses = true;
        return OperationStatus::success;
      }
    }
  }
  return OperationStatus::success;
}

inline OperationStatus generateChildrenFromSource(const Polymap* pMap,
                                           int owner_index,
                                           const Point2d& source,
                                           const exact_geometry::Point& exact_source,
                                           double start_angle,
                                           double end_angle,
                                           double gcost,
                                           const VisibilityRegion& visibility_region,
                                           std::vector<Child>& output,
                                           const StopToken& stop_token,
                                           std::string* error) {
  output.clear();
  if (error)
    error->clear();
  const auto fail = [&](const std::string& message) {
    output.clear();
    if (error)
      *error = message;
    return OperationStatus::failure;
  };
  const auto stopped = [&]() {
    output.clear();
    if (error)
      error->clear();
    return OperationStatus::stopped;
  };
  if (stop_token.poll())
    return stopped();
  if (!pMap)
    return fail("Cannot generate children without a Polymap");
  if (visibility_region.size() < 2)
    return fail("Visibility region has fewer than two vertices");

  std::string validation_error;
  const auto visibility_status =
    pMap->validateVisibilityRegion(visibility_region, stop_token, &validation_error);
  if (visibility_status == OperationStatus::stopped)
    return stopped();
  if (visibility_status == OperationStatus::failure)
    return fail(validation_error);

  std::vector<double> theta_list(visibility_region.size());
  for (size_t i = 0; i < visibility_region.size(); ++i) {
    if (stop_token.poll())
      return stopped();
    const auto& position = visibility_region[i].position;
    theta_list[i] = atan2(position.second - source.second, position.first - source.first);
    theta_list[i] = start_angle + normalize_angle_positive(theta_list[i] - start_angle);
  }

  std::vector<bool> is_a_left_gap(visibility_region.size(), false);
  std::vector<size_t> valid_child_indices;
  for (size_t i = 0; i + 1 < visibility_region.size(); ++i) {
    if (stop_token.poll())
      return stopped();
    const size_t next = i + 1;
    const auto& current_endpoint = visibility_region[i];
    const auto& next_endpoint = visibility_region[next];
    if (exactPoint(current_endpoint) == exactPoint(next_endpoint))
      continue;
    if (exact_geometry::isSameDirectedRay(
          exact_source, exactPoint(current_endpoint), exactPoint(next_endpoint))) {
      if (std::holds_alternative<std::monostate>(current_endpoint.support) &&
          std::holds_alternative<std::monostate>(next_endpoint.support)) {
        continue;
      }
      bool supports_consecutive = false;
      const auto support_status = pMap->boundarySupportsConsecutive(
        next_endpoint, current_endpoint, supports_consecutive, stop_token);
      if (support_status == OperationStatus::stopped)
        return stopped();
      if (support_status == OperationStatus::failure)
        return fail("Visibility boundary support validation failed");
      if (!supports_consecutive) {
        const auto current_distance =
          (exactPoint(current_endpoint) - exact_source).squared_length();
        const auto next_distance = (exactPoint(next_endpoint) - exact_source).squared_length();
        if (current_distance > next_distance)
          is_a_left_gap[i] = true;
        valid_child_indices.emplace_back(i);
      }
    }
  }

  struct PendingChild {
    BoundaryEndpoint near_endpoint;
    BoundaryEndpoint far_endpoint;
    bool is_left_gap = false;
    double start_angle = 0.0;
    double end_angle = 0.0;
  };
  std::vector<PendingChild> pending_children;
  pending_children.reserve(valid_child_indices.size());
  for (const auto i : valid_child_indices) {
    if (stop_token.poll())
      return stopped();
    const size_t next = i + 1;
    const size_t near_index = is_a_left_gap[i] ? next : i;
    const size_t far_index = is_a_left_gap[i] ? i : next;
    const auto& near_endpoint = visibility_region[near_index];
    const auto& far_endpoint = visibility_region[far_index];
    const auto* near_vertex = std::get_if<ObstacleVertexId>(&near_endpoint.support);
    if (!near_vertex) {
      return fail("Visibility gap has an edge-interior or unsupported near anchor at position " +
                  std::to_string(near_index));
    }
    const std::pair<int, int> near_topology = {near_vertex->obstacle, near_vertex->vertex};
    if (!pMap->isValidTopology(near_topology))
      return fail("Visibility gap has an invalid near anchor topology");

    double child_start_angle = 0.0;
    double child_end_angle = 0.0;
    if (is_a_left_gap[i]) {
      const auto next_obs = pMap->getNextObs(near_topology);
      if (!next_obs)
        return fail("Cannot find the next obstacle vertex for a visibility gap");
      child_start_angle = normalize_angle(theta_list[next]);
      const double contour_angle = atan2(next_obs->second - near_endpoint.position.second,
                                         next_obs->first - near_endpoint.position.first);
      child_end_angle =
        child_start_angle + normalize_angle_positive(contour_angle - child_start_angle);
    } else {
      const auto prev_obs = pMap->getPrevObs(near_topology);
      if (!prev_obs)
        return fail("Cannot find the previous obstacle vertex for a visibility gap");
      const double contour_angle = atan2(prev_obs->second - near_endpoint.position.second,
                                         prev_obs->first - near_endpoint.position.first);
      child_start_angle = contour_angle;
      child_end_angle = contour_angle + normalize_angle_positive(theta_list[i] - contour_angle);
    }

    // The near endpoint is topologically an obstacle vertex.  Recover the
    // integer coordinate from that identity instead of rounding a projected
    // double, so the Child contract remains strictly integral even for a
    // continuous root source.
    const auto& canonical = pMap->obstacles()[static_cast<size_t>(near_topology.first)]
                              .ordered_vertices_[static_cast<size_t>(near_topology.second)];
    BoundaryEndpoint canonical_near_endpoint = near_endpoint;
    canonical_near_endpoint.position = {static_cast<double>(canonical.first),
                                        static_cast<double>(canonical.second)};
    canonical_near_endpoint.exact_position =
      exact_geometry::Point(canonical.first, canonical.second);
    pending_children.push_back(PendingChild{std::move(canonical_near_endpoint),
                                            far_endpoint,
                                            is_a_left_gap[i],
                                            child_start_angle,
                                            child_end_angle});
  }

  std::vector<Child> generated_children;
  generated_children.reserve(pending_children.size());
  for (size_t i = 0; i < pending_children.size(); ++i) {
    if (stop_token.poll())
      return stopped();
    auto& pending = pending_children[i];
    const double child_gcost =
      gcost + std::hypot(source.first - pending.near_endpoint.position.first,
                         source.second - pending.near_endpoint.position.second);
    generated_children.emplace_back(owner_index,
                                    static_cast<int>(i),
                                    pending.near_endpoint,
                                    pending.far_endpoint,
                                    pending.is_left_gap,
                                    pending.start_angle,
                                    pending.end_angle,
                                    child_gcost);
  }
  if (stop_token.poll())
    return stopped();
  output = std::move(generated_children);
  (void)end_angle;  // retained in the helper signature for the Node contract
  return OperationStatus::success;
}


}  // namespace core_impl

}  // namespace raystar

#endif  // RAYSTAR_CORE_DETAIL_H_
