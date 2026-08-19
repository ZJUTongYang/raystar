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

RaystarNode::RaystarNode(const rclcpp::NodeOptions& options) : Node("raystar", options) {
  // Register before declaration so command-line/YAML overrides are checked by
  // the same validator as runtime updates.  The callback has no side effects:
  // rclcpp owns the all-or-none parameter commit after every callback accepts.
  parameter_validation_callback_ = add_on_set_parameters_callback(validateParameterChanges);
  for (const auto& spec : kIntegerParameterSpecs) {
    declare_parameter<int64_t>(spec.name, spec.default_value, makeParameterDescriptor(spec));
  }

  rcl_interfaces::msg::ParameterDescriptor map_topic_descriptor;
  map_topic_descriptor.description =
    "OccupancyGrid topic cached by the Action server. This startup-only "
    "parameter must match the map observed by Action clients.";
  map_topic_descriptor.read_only = true;
  declare_parameter<std::string>("map_topic", "/map", map_topic_descriptor);

  rcl_interfaces::msg::ParameterDescriptor legacy_service_descriptor;
  legacy_service_descriptor.description =
    "Expose the compatibility Service whose request still embeds a complete "
    "OccupancyGrid. Prefer the cached-map Action for bounded DDS traffic.";
  legacy_service_descriptor.read_only = true;
  declare_parameter<bool>("enable_legacy_map_service", true, legacy_service_descriptor);

  RequestConfiguration configuration;
  std::string configuration_error;
  if (!loadRequestConfiguration(configuration, configuration_error)) {
    throw rclcpp::exceptions::InvalidParameterValueException(
      "Invalid Raystar server configuration: " + configuration_error);
  }

  if (get_parameter("enable_legacy_map_service").as_bool()) {
    service_ = create_service<raystar_interfaces::srv::GetRaystarPaths>(
      "~/get_raystar_paths",
      std::bind(&RaystarNode::handleService, this, std::placeholders::_1, std::placeholders::_2));
    RCLCPP_WARN(get_logger(),
                "Legacy full-map Service is enabled; prefer the cached-map Action and "
                "disable enable_legacy_map_service in resource-constrained deployments");
  }

  map_status_publisher_ = create_publisher<raystar_interfaces::msg::MapStatus>(
    "~/map_status", rclcpp::QoS(rclcpp::KeepLast(1)).transient_local());
  const std::string map_topic = get_parameter("map_topic").as_string();
  map_subscription_ = create_subscription<nav_msgs::msg::OccupancyGrid>(
    map_topic,
    rclcpp::QoS(rclcpp::KeepLast(1)).transient_local().reliable(),
    [this](nav_msgs::msg::OccupancyGrid::ConstSharedPtr map) { handleMap(std::move(map)); });

  action_server_ = rclcpp_action::create_server<PlanAction>(
    this,
    "~/plan_paths",
    std::bind(&RaystarNode::handleActionGoal, this, std::placeholders::_1, std::placeholders::_2),
    std::bind(&RaystarNode::handleActionCancel, this, std::placeholders::_1),
    std::bind(&RaystarNode::handleActionAccepted, this, std::placeholders::_1));
  goal_set_action_server_ = rclcpp_action::create_server<GoalSetAction>(
    this,
    "~/plan_goal_set",
    std::bind(
      &RaystarNode::handleGoalSetActionGoal, this, std::placeholders::_1, std::placeholders::_2),
    std::bind(&RaystarNode::handleGoalSetActionCancel, this, std::placeholders::_1),
    std::bind(&RaystarNode::handleGoalSetActionAccepted, this, std::placeholders::_1));
  transition_action_server_ = rclcpp_action::create_server<TransitionAction>(
    this,
    "~/build_transition_graph",
    std::bind(
      &RaystarNode::handleTransitionActionGoal, this, std::placeholders::_1, std::placeholders::_2),
    std::bind(&RaystarNode::handleTransitionActionCancel, this, std::placeholders::_1),
    std::bind(&RaystarNode::handleTransitionActionAccepted, this, std::placeholders::_1));

  auto latched = rclcpp::QoS(rclcpp::KeepLast(1)).transient_local();
  non_homotopic_pub_ =
    create_publisher<visualization_msgs::msg::MarkerArray>("~/non_homotopic_paths", latched);
  poly_obstacle_pub_ =
    create_publisher<visualization_msgs::msg::MarkerArray>("~/poly_obstacles", latched);
  debug_tree_pub_ = create_publisher<visualization_msgs::msg::MarkerArray>("~/debug_tree", latched);
  cdt_pub_ = create_publisher<visualization_msgs::msg::MarkerArray>("~/cdt", latched);

  const int64_t path_republish_period_ms =
    get_parameter("path_visualization_republish_period_ms").as_int();
  if (path_republish_period_ms > 0) {
    path_visualization_timer_ =
      create_wall_timer(std::chrono::milliseconds(path_republish_period_ms),
                        std::bind(&RaystarNode::republishCachedPathVisualization, this));
  }

  RCLCPP_INFO(get_logger(),
              "Raystar node initialized: occupied_threshold=%d "
              "max_k=%d max_cost_bounded_paths=%zu max_multi_goal_count=%zu "
              "max_transition_configurations=%zu max_transition_pairs=%zu max_nodes=%zu "
              "planning_timeout_ms=%lld max_map_cells=%zu max_map_bytes=%zu "
              "max_path_points=%zu max_debug_nodes=%zu max_response_bytes=%zu "
              "path_visualization_republish_period_ms=%lld map_topic=%s",
              configuration.occupied_threshold,
              configuration.planning.max_k,
              configuration.planning.max_cost_bounded_paths,
              configuration.planning.max_multi_goal_count,
              configuration.planning.max_transition_configurations,
              configuration.planning.max_transition_pairs,
              configuration.planning.max_nodes,
              static_cast<long long>(configuration.planning.planning_timeout.count()),
              configuration.planning.max_map_cells,
              configuration.planning.max_map_bytes,
              configuration.planning.max_path_points,
              configuration.planning.max_debug_nodes,
              configuration.planning.max_response_bytes,
              static_cast<long long>(path_republish_period_ms),
              map_topic.c_str());

  // Start the persistent capacity-one Action worker only after every other
  // constructor operation that can throw has completed.  Accepted/cancel
  // callbacks only update bounded state and notify this worker.
  action_worker_ = std::thread(&RaystarNode::actionWorkerLoop, this);
}

RaystarNode::~RaystarNode() {
  if (path_visualization_timer_)
    path_visualization_timer_->cancel();
  shutting_down_.store(true, std::memory_order_release);
  {
    std::lock_guard<std::mutex> lock(action_worker_mutex_);
    stop_action_worker_ = true;
  }
  action_worker_cv_.notify_one();
  if (action_worker_.joinable())
    action_worker_.join();

  // main() destroys the node before rclcpp::shutdown(), so publishers are
  // still usable during a normal SIGINT/clean exit. Serialize this final
  // best-effort DELETEALL against any last timer callback and release retained
  // search/Marker state. clearVisualizationsLocked() contains its own
  // per-publisher exception barriers for middleware teardown races.
  std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
  clearVisualizationsLocked();
}

bool RaystarNode::loadRequestConfiguration(RequestConfiguration& configuration,
                                           std::string& error) const {
  const auto parameters = get_parameters({"max_k",
                                          "max_cost_bounded_paths",
                                          "max_multi_goal_count",
                                          "max_transition_configurations",
                                          "max_transition_pairs",
                                          "max_nodes",
                                          "planning_timeout_ms",
                                          "max_map_cells",
                                          "max_map_bytes",
                                          "max_path_points",
                                          "max_debug_nodes",
                                          "max_response_bytes",
                                          "occupied_threshold"});
  if (parameters.size() != 13) {
    error = "could not read the complete ROS server parameter set";
    return false;
  }
  const auto validation = validateParameterChanges(parameters);
  if (!validation.successful) {
    error = validation.reason;
    return false;
  }

  const int64_t max_k = parameters[0].as_int();
  const int64_t max_cost_bounded_paths = parameters[1].as_int();
  const int64_t max_multi_goal_count = parameters[2].as_int();
  const int64_t max_transition_configurations = parameters[3].as_int();
  const int64_t max_transition_pairs = parameters[4].as_int();
  const int64_t max_nodes = parameters[5].as_int();
  const int64_t planning_timeout_ms = parameters[6].as_int();
  const int64_t max_map_cells = parameters[7].as_int();
  const int64_t max_map_bytes = parameters[8].as_int();
  const int64_t max_path_points = parameters[9].as_int();
  const int64_t max_debug_nodes = parameters[10].as_int();
  const int64_t max_response_bytes = parameters[11].as_int();
  const int64_t occupied_threshold = parameters[12].as_int();
  configuration = RequestConfiguration{};
  configuration.occupied_threshold = static_cast<int>(occupied_threshold);
  configuration.planning.max_k = static_cast<int>(max_k);
  configuration.planning.max_cost_bounded_paths = static_cast<size_t>(max_cost_bounded_paths);
  configuration.planning.max_multi_goal_count = static_cast<size_t>(max_multi_goal_count);
  configuration.planning.max_transition_configurations =
    static_cast<size_t>(max_transition_configurations);
  configuration.planning.max_transition_pairs = static_cast<size_t>(max_transition_pairs);
  configuration.planning.max_nodes = static_cast<size_t>(max_nodes);
  configuration.planning.planning_timeout = std::chrono::milliseconds(planning_timeout_ms);
  configuration.planning.max_map_cells = static_cast<size_t>(max_map_cells);
  configuration.planning.max_map_bytes = static_cast<size_t>(max_map_bytes);
  configuration.planning.max_path_points = static_cast<size_t>(max_path_points);
  configuration.planning.max_debug_nodes = static_cast<size_t>(max_debug_nodes);
  configuration.planning.max_response_bytes = static_cast<size_t>(max_response_bytes);
  error.clear();
  return true;
}

void RaystarNode::handleMap(nav_msgs::msg::OccupancyGrid::ConstSharedPtr map) {
  if (!map || shutting_down_.load(std::memory_order_acquire))
    return;

  try {
    RequestConfiguration configuration;
    std::string error;
    if (!loadRequestConfiguration(configuration, error)) {
      RCLCPP_ERROR(get_logger(),
                   "Cannot admit cached map because server configuration is invalid: %s",
                   error.c_str());
      return;
    }

    // Run the same complete validation and allocation budget used by planning
    // before making this DDS-owned immutable sample visible to Action goals.
    // The temporary binary map is released immediately; request execution
    // converts the retained raw occupancy values using its own allow_unknown
    // snapshot.
    GridMap validation_map;
    const StopPredicate stop_requested = [this]() {
      return shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok();
    };
    if (!occupancyGridToBinaryMap(
          *map, false, configuration, stop_requested, validation_map, error)) {
      RCLCPP_WARN(get_logger(), "Rejected map cache update: %s", error.c_str());
      return;
    }

    const raystar_interfaces::MapId id = raystar_interfaces::computeMapId(*map);
    std::uint64_t generation = 0;
    {
      std::lock_guard<std::mutex> lock(map_cache_mutex_);
      generation = cached_map_.generation == std::numeric_limits<std::uint64_t>::max()
                     ? 1U
                     : cached_map_.generation + 1U;
      cached_map_.message = map;
      cached_map_.id = id;
      cached_map_.generation = generation;
      // Even an identical republished OccupancyGrid is a new admitted
      // snapshot. Drop the one strong environment entry while holding the map
      // lock so an older in-flight GCP cannot race this update and repopulate
      // stale geometry afterwards.
      transition_environment_cache_.clear();
    }
    raystar_interfaces::msg::MapStatus status;
    status.header = map->header;
    status.map_id = id;
    status.generation = generation;
    status.width = map->info.width;
    status.height = map->info.height;
    status.resolution = map->info.resolution;
    status.occupied_threshold = static_cast<uint32_t>(configuration.occupied_threshold);
    status.environment_identity_version = raystar_interfaces::kEnvironmentIdentityVersion;
    status.occupancy_semantics_version = raystar_interfaces::kOccupancySemanticsVersion;
    status.geometry_semantics_version = raystar_interfaces::kGeometrySemanticsVersion;
    status.topology_semantics_version = raystar_interfaces::kTopologySemanticsVersion;
    status.environment_id_disallow_unknown = raystar_interfaces::computeEnvironmentId(
      *map, configuration.occupied_threshold, false);
    status.environment_id_allow_unknown = raystar_interfaces::computeEnvironmentId(
      *map, configuration.occupied_threshold, true);
    map_status_publisher_->publish(status);
    RCLCPP_INFO(get_logger(),
                "Cached admitted OccupancyGrid generation=%llu",
                static_cast<unsigned long long>(generation));
  } catch (const std::exception& exception) {
    RCLCPP_ERROR(get_logger(), "Could not cache OccupancyGrid: %s", exception.what());
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Could not cache OccupancyGrid because of an unknown exception");
  }
}

bool RaystarNode::resolveCachedMap(const raystar_interfaces::MapId& requested_id,
                                   nav_msgs::msg::OccupancyGrid::ConstSharedPtr& map,
                                   std::string& error) const {
  map.reset();
  if (raystar_interfaces::isZeroMapId(requested_id)) {
    error = "Invalid map_id: cached-map Action requires a non-zero map identity";
    return false;
  }

  std::lock_guard<std::mutex> lock(map_cache_mutex_);
  if (!cached_map_.message) {
    error = "No admitted OccupancyGrid is cached; publish the configured map_topic first";
    return false;
  }
  if (!raystar_interfaces::mapIdsEqual(requested_id, cached_map_.id)) {
    error = "Requested map_id is not the current cached OccupancyGrid; refresh the map snapshot";
    return false;
  }
  map = cached_map_.message;
  error.clear();
  return true;
}

std::shared_ptr<const Polymap> RaystarNode::findCachedTransitionEnvironment(
  const nav_msgs::msg::OccupancyGrid& grid,
  const raystar_interfaces::MapId& map_id,
  bool allow_unknown,
  const RequestConfiguration& configuration,
  const PolymapEndpoint& base,
  const std::vector<PolymapEndpoint>& goals) {
  std::lock_guard<std::mutex> lock(map_cache_mutex_);
  if (!cached_map_.message || cached_map_.message.get() != &grid ||
      !raystar_interfaces::mapIdsEqual(cached_map_.id, map_id)) {
    return {};
  }
  const TransitionEnvironmentPolicy policy{allow_unknown,
                                           configuration.occupied_threshold,
                                           configuration.planning.max_map_cells,
                                           configuration.planning.max_map_bytes};
  std::vector<TransitionEnvironmentEndpoint> cache_goals;
  cache_goals.reserve(goals.size());
  for (const auto& goal : goals) cache_goals.emplace_back(goal);
  return transition_environment_cache_.find(TransitionEnvironmentKey(
    cached_map_.generation, policy, TransitionEnvironmentEndpoint(base), std::move(cache_goals)));
}

void RaystarNode::cacheCompletedTransitionEnvironment(const nav_msgs::msg::OccupancyGrid& grid,
                                                      const raystar_interfaces::MapId& map_id,
                                                      bool allow_unknown,
                                                      const RequestConfiguration& configuration,
                                                      const PolymapEndpoint& base,
                                                      const std::vector<PolymapEndpoint>& goals,
                                                      std::shared_ptr<const Polymap> environment) {
  if (!environment)
    return;
  std::lock_guard<std::mutex> lock(map_cache_mutex_);
  if (!cached_map_.message || cached_map_.message.get() != &grid ||
      !raystar_interfaces::mapIdsEqual(cached_map_.id, map_id)) {
    return;
  }
  const TransitionEnvironmentPolicy policy{allow_unknown,
                                           configuration.occupied_threshold,
                                           configuration.planning.max_map_cells,
                                           configuration.planning.max_map_bytes};
  std::vector<TransitionEnvironmentEndpoint> cache_goals;
  cache_goals.reserve(goals.size());
  for (const auto& goal : goals) cache_goals.emplace_back(goal);
  transition_environment_cache_.store(
    TransitionEnvironmentKey(
      cached_map_.generation, policy, TransitionEnvironmentEndpoint(base), std::move(cache_goals)),
    std::move(environment));
}

bool RaystarNode::occupancyGridToBinaryMap(const nav_msgs::msg::OccupancyGrid& grid,
                                           bool allow_unknown,
                                           const RequestConfiguration& configuration,
                                           const StopPredicate& stop_requested,
                                           GridMap& output,
                                           std::string& error) const {
  output = GridMap{};
  error.clear();
  const size_t width = static_cast<size_t>(grid.info.width);
  const size_t height = static_cast<size_t>(grid.info.height);

  if (grid.header.frame_id.empty()) {
    error = "Invalid map: header.frame_id must not be empty";
    return false;
  }

  if (!std::isfinite(static_cast<double>(grid.info.resolution)) || grid.info.resolution <= 0.0f) {
    error = "Invalid map: resolution must be finite and greater than zero";
    return false;
  }

  const auto& origin = grid.info.origin;
  if (!std::isfinite(origin.position.x) || !std::isfinite(origin.position.y) ||
      !std::isfinite(origin.position.z)) {
    error = "Invalid map: origin position must contain only finite coordinates";
    return false;
  }
  if (origin.position.z != 0.0) {
    error = "Invalid map: origin z coordinate must be zero for 2D planning";
    return false;
  }

  const auto& orientation = origin.orientation;
  if (!std::isfinite(orientation.x) || !std::isfinite(orientation.y) ||
      !std::isfinite(orientation.z) || !std::isfinite(orientation.w)) {
    error = "Invalid map: origin orientation must contain only finite values";
    return false;
  }
  const double quaternion_norm_squared =
    orientation.x * orientation.x + orientation.y * orientation.y + orientation.z * orientation.z +
    orientation.w * orientation.w;
  if (!std::isfinite(quaternion_norm_squared) ||
      std::abs(quaternion_norm_squared - 1.0) > kQuaternionNormTolerance) {
    error = "Invalid map: origin orientation must be a unit quaternion";
    return false;
  }
  if (orientation.x != 0.0 || orientation.y != 0.0 || orientation.z != 0.0) {
    error =
      "Invalid map: origin orientation must be the identity rotation; "
      "rotated or tilted occupancy grids are not supported";
    return false;
  }

  // Reuse Core's allocation-free admission check so ROS and direct callers
  // cannot drift on overflow, signed-indexing, or byte-budget semantics.
  MapResourceEstimate map_estimate;
  if (!validateMapResourceBudget(
        width, height, grid.data.size(), configuration.planning, map_estimate, error)) {
    return false;
  }
  const size_t cell_count = map_estimate.cell_count;

  const double metric_resolution = static_cast<double>(grid.info.resolution);
  const double extent_x = gridCostToMetric(static_cast<double>(width), grid.info.resolution);
  const double extent_y = gridCostToMetric(static_cast<double>(height), grid.info.resolution);
  const double max_world_x = std::fma(1.0, extent_x, origin.position.x);
  const double max_world_y = std::fma(1.0, extent_y, origin.position.y);
  if (!std::isfinite(max_world_x) || !std::isfinite(max_world_y)) {
    error = "Invalid map: world-coordinate extent is not finite";
    return false;
  }
  if (!hasMetricFaithfulWorldTransform(
        origin.position.x, origin.position.y, extent_x, extent_y, metric_resolution)) {
    error =
      "Invalid map: binary64 grid-to-world transform lacks metric fidelity "
      "for published path lengths";
    return false;
  }

  GridMap map;
  map.width = grid.info.width;
  map.height = grid.info.height;
  map.resolution = grid.info.resolution;
  map.origin_x = grid.info.origin.position.x;
  map.origin_y = grid.info.origin.position.y;

  map.data.resize(cell_count);
  for (size_t i = 0; i < cell_count; ++i) {
    if ((i & 0xfffu) == 0u && stop_requested && stop_requested()) {
      error = "Planning was cancelled while converting the occupancy grid";
      return false;
    }
    const int value = static_cast<int>(grid.data[i]);
    if (value < -1 || value > 100) {
      const size_t x = i % width;
      const size_t y = i / width;
      error = "Invalid map: occupancy value " + std::to_string(value) + " at cell (" +
              std::to_string(x) + ", " + std::to_string(y) +
              ") is outside the supported range [-1, 100]";
      return false;
    }
    if (value == -1)
      map.data[i] = allow_unknown ? 0 : 1;
    else
      map.data[i] = value >= configuration.occupied_threshold ? 1 : 0;
  }
  output = std::move(map);
  return true;
}

void RaystarNode::handleService(
  const std::shared_ptr<raystar_interfaces::srv::GetRaystarPaths::Request> request,
  std::shared_ptr<raystar_interfaces::srv::GetRaystarPaths::Response> response) {
  resetPlanningResponse(*response);
  setRequestedPathCount(*response, request ? request->k : 0);
  if (request) {
    response->result_info.search_mode = request->search_mode;
    response->result_info.requested_max_path_length = request->max_path_length;
    response->result_info.debug_requested = request->include_debug;
    response->result_info.debug_output_complete = !request->include_debug;
  }
  if (!request) {
    response->result_info.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
    response->message = "Raystar Service request is unavailable";
    return;
  }
  bool expected_idle = false;
  if (!planning_busy_.compare_exchange_strong(expected_idle, true, std::memory_order_acq_rel)) {
    response->result_info.status = PlanningResultInfo::STATUS_BUSY;
    response->message = "Raystar planner is busy with another Service request or Action goal";
    return;
  }

  PlanningSlotRelease release_slot(planning_busy_);
  const auto map_id = raystar_interfaces::computeMapId(request->map);
  executePlanning(*request, *response, request->map, map_id, [this]() {
    return shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok();
  });
}

rclcpp_action::GoalResponse RaystarNode::handleActionGoal(
  const rclcpp_action::GoalUUID& uuid, std::shared_ptr<const PlanAction::Goal> goal) {
  if (!goal || shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok()) {
    return rclcpp_action::GoalResponse::REJECT;
  }

  std::shared_ptr<std::atomic<bool>> cancel_requested;
  try {
    cancel_requested = std::make_shared<std::atomic<bool>>(false);
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Rejecting Action goal: could not allocate cancellation state");
    return rclcpp_action::GoalResponse::REJECT;
  }

  bool expected_idle = false;
  if (!planning_busy_.compare_exchange_strong(expected_idle, true, std::memory_order_acq_rel)) {
    RCLCPP_WARN(get_logger(), "Rejecting Action goal because the capacity-one planner is busy");
    return rclcpp_action::GoalResponse::REJECT;
  }

  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    active_goal_id_ = uuid;
    active_goal_reserved_ = true;
    active_goal_cancel_ = std::move(cancel_requested);
  }
  return rclcpp_action::GoalResponse::ACCEPT_AND_EXECUTE;
}

rclcpp_action::CancelResponse RaystarNode::handleActionCancel(
  const std::shared_ptr<PlanGoalHandle> goal_handle) {
  if (!goal_handle)
    return rclcpp_action::CancelResponse::REJECT;

  std::lock_guard<std::mutex> lock(action_state_mutex_);
  if (!active_goal_reserved_ || active_goal_id_ != goal_handle->get_goal_id() ||
      !active_goal_cancel_) {
    return rclcpp_action::CancelResponse::REJECT;
  }
  // Core polls this per-goal flag directly.  rclcpp_action transitions the
  // goal to CANCELING immediately after this callback returns ACCEPT; the
  // worker waits for that transition before publishing its canceled result.
  active_goal_cancel_->store(true, std::memory_order_release);
  return rclcpp_action::CancelResponse::ACCEPT;
}

}  // namespace raystar
