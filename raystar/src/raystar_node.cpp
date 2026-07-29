#include "raystar_node.h"

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
#include <iomanip>
#include <utility>
#include <exception>

namespace raystar
{

namespace
{

// This tolerance validates quaternion encoding only. Rotation components are
// still required to be exactly zero, so no geometric yaw/tilt is ignored.
constexpr double kQuaternionNormTolerance = 1e-12;
constexpr int64_t kDefaultOccupiedThreshold = 99;
constexpr int64_t kDefaultMaxK = 100;
constexpr int64_t kDefaultMaxNodes = 10000;
constexpr int64_t kDefaultPlanningTimeoutMs = 5000;
constexpr int64_t kDefaultMaxMapCells =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxMapCells);
constexpr int64_t kDefaultMaxMapBytes =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxMapBytes);
constexpr int64_t kDefaultMaxPathPoints =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxPathPoints);
constexpr int64_t kDefaultMaxDebugNodes =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxDebugNodes);
constexpr int64_t kDefaultMaxResponseBytes =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxResponseBytes);
constexpr int64_t kDefaultPathVisualizationRepublishPeriodMs = 2000;
constexpr int64_t kMaxTimerPeriodMs =
  std::chrono::duration_cast<std::chrono::milliseconds>(
    std::chrono::nanoseconds::max() - std::chrono::milliseconds(1)).count();
constexpr int64_t kMaxIntParameterValue =
  static_cast<int64_t>(std::numeric_limits<int>::max());
constexpr int64_t kMaxPlanningTimeoutMs =
  std::chrono::milliseconds::max().count() - 1;

constexpr int64_t maxSizeTParameterValue()
{
  if constexpr (sizeof(size_t) < sizeof(int64_t))
    return static_cast<int64_t>(std::numeric_limits<size_t>::max());
  return std::numeric_limits<int64_t>::max();
}

// Reserve a small fixed amount for CDR sequence/string headers.  The message
// field is charged separately (including every bounded notice appended below)
// so a small, valid max_response_bytes value does not get consumed twice by a
// diagnostic reserve.  Per-element estimates deliberately include alignment
// and headroom, rather than pretending to be an exact RMW serialization size.
constexpr size_t kMinimumResponseBytes = 1024;
constexpr size_t kResponseBaseBytes = 256;
constexpr size_t kPathBaseBytes = 128;
constexpr size_t kPoseBaseBytes = 96;
constexpr size_t kDebugNodeBytes = 32;
constexpr size_t kMaxFrameIdBytes = 256;
constexpr size_t kMaxDiagnosticBytes = 768;

struct IntegerParameterSpec
{
  const char* name;
  int64_t default_value;
  int64_t minimum;
  int64_t maximum;
  const char* description;
  bool read_only;
};

constexpr std::array<IntegerParameterSpec, 10> kIntegerParameterSpecs{{
  {"occupied_threshold", kDefaultOccupiedThreshold, 1, 100,
    "Known OccupancyGrid values at or above this threshold are occupied.",
    false},
  {"max_k", kDefaultMaxK, 1, kMaxIntParameterValue,
    "Largest requested number of paths accepted by the server.", false},
  {"max_nodes", kDefaultMaxNodes, 1, kMaxIntParameterValue,
    "Maximum fully constructed search nodes, including the root.", false},
  {"planning_timeout_ms", kDefaultPlanningTimeoutMs, 1,
    kMaxPlanningTimeoutMs,
    "Cooperative planning deadline in milliseconds.", false},
  {"max_map_cells", kDefaultMaxMapCells, 1, kMaxIntParameterValue,
    "Maximum OccupancyGrid width multiplied by height.", false},
  {"max_map_bytes", kDefaultMaxMapBytes,
    static_cast<int64_t>(kEstimatedPlannerMapBytesPerCell),
    maxSizeTParameterValue(),
    "Conservative map working-set admission budget in bytes.", false},
  {"max_path_points", kDefaultMaxPathPoints, 2, kMaxIntParameterValue,
    "Maximum total sampled path poses in one response.", false},
  {"max_debug_nodes", kDefaultMaxDebugNodes, 0, kMaxIntParameterValue,
    "Maximum structured debug nodes; zero disables debug-node output.", false},
  {"max_response_bytes", kDefaultMaxResponseBytes,
    static_cast<int64_t>(kMinimumResponseBytes), maxSizeTParameterValue(),
    "Conservative response and visualization payload budget in bytes.", false},
  {"path_visualization_republish_period_ms",
    kDefaultPathVisualizationRepublishPeriodMs, 0, kMaxTimerPeriodMs,
    "Startup-only cached path MarkerArray republish period; zero disables it.",
    true},
}};

const IntegerParameterSpec* findIntegerParameterSpec(const std::string& name)
{
  const auto found = std::find_if(
    kIntegerParameterSpecs.begin(), kIntegerParameterSpecs.end(),
    [&name](const IntegerParameterSpec& spec) {
      return name == spec.name;
    });
  return found == kIntegerParameterSpecs.end() ? nullptr : &*found;
}

rcl_interfaces::msg::ParameterDescriptor makeParameterDescriptor(
  const IntegerParameterSpec& spec)
{
  rcl_interfaces::msg::IntegerRange range;
  range.from_value = spec.minimum;
  range.to_value = spec.maximum;
  range.step = 1;

  rcl_interfaces::msg::ParameterDescriptor descriptor;
  descriptor.description = spec.description;
  descriptor.read_only = spec.read_only;
  descriptor.integer_range.emplace_back(range);
  descriptor.additional_constraints = spec.read_only ?
    "May be overridden at startup but cannot be changed while the node runs." :
    "An accepted update affects only requests that take a later configuration snapshot.";
  return descriptor;
}

rcl_interfaces::msg::SetParametersResult validateParameterChanges(
  const std::vector<rclcpp::Parameter>& parameters)
{
  rcl_interfaces::msg::SetParametersResult result;
  result.successful = true;
  for (const auto& parameter : parameters)
  {
    const IntegerParameterSpec* spec =
      findIntegerParameterSpec(parameter.get_name());
    // Other components may register their own parameters on the node.  Do not
    // reject values outside Raystar's owned parameter set.
    if (spec == nullptr)
      continue;
    if (parameter.get_type() != rclcpp::ParameterType::PARAMETER_INTEGER)
    {
      result.successful = false;
      result.reason = "parameter '" + parameter.get_name() +
        "' must be an integer";
      return result;
    }
    const int64_t value = parameter.as_int();
    if (value < spec->minimum || value > spec->maximum)
    {
      result.successful = false;
      result.reason = "parameter '" + parameter.get_name() +
        "' must be between " + std::to_string(spec->minimum) + " and " +
        std::to_string(spec->maximum);
      return result;
    }
  }
  return result;
}

using PlanningResultInfo = raystar_interfaces::msg::PlanningResultInfo;

// MarkerArray topics are separate from the planning result, but they are
// still built synchronously in the admitted planning execution. Give each of the four
// topics a quarter of the configured response budget and use conservative
// per-element estimates before growing ROS message vectors.  This is an
// admission guard for serialization allocations, not an exact DDS wire-size
// calculation.
constexpr size_t kVisualizationTopicCount = 4;
constexpr size_t kMarkerPointBytes = 64;
constexpr size_t kMarkerEntryBytes = 512;

size_t markerTopicBudget(size_t total_budget)
{
  return std::max<size_t>(1, total_budget / kVisualizationTopicCount);
}

size_t markerPointLimit(size_t topic_budget)
{
  return topic_budget / kMarkerPointBytes;
}

size_t markerEntryLimit(size_t topic_budget)
{
  return topic_budget / kMarkerEntryBytes;
}

bool canAppendCount(size_t current, size_t addition, size_t limit)
{
  return addition <= limit && current <= limit - addition;
}

visualization_msgs::msg::MarkerArray makeMarkerSnapshot(
  const std::string& frame_id, const rclcpp::Time& stamp)
{
  visualization_msgs::msg::MarkerArray array;
  visualization_msgs::msg::Marker clear_marker;
  clear_marker.header.frame_id = frame_id;
  clear_marker.header.stamp = stamp;
  clear_marker.action = visualization_msgs::msg::Marker::DELETEALL;
  array.markers.emplace_back(std::move(clear_marker));
  return array;
}

bool validatePlanarPose(
  const geometry_msgs::msg::PoseStamped& pose, const std::string& label,
  const std::string& map_frame, std::string& error)
{
  if (pose.header.frame_id.empty())
  {
    error = label + " frame_id must not be empty";
    return false;
  }
  if (pose.header.frame_id != map_frame)
  {
    error = label + " frame_id '" + pose.header.frame_id +
      "' does not match map frame_id '" + map_frame + "'";
    return false;
  }

  const auto& position = pose.pose.position;
  if (!std::isfinite(position.x) || !std::isfinite(position.y) ||
      !std::isfinite(position.z))
  {
    error = label + " position must contain only finite coordinates";
    return false;
  }
  if (position.z != 0.0)
  {
    error = label + " z coordinate must be zero for 2D planning";
    return false;
  }
  return true;
}

// Path waypoints are expressed in the continuous grid coordinate system.
// Unlike coordinate_utils::continuousMapToWorld(), this conversion also
// accepts obstacle vertices on the outer map boundary (x == width or
// y == height), which are valid geometric waypoints even though they are not
// cell interiors.
bool continuousGridToWorld(
  const GridMap& grid_map, const Point2d& point, double& wx, double& wy)
{
  wx = std::numeric_limits<double>::quiet_NaN();
  wy = std::numeric_limits<double>::quiet_NaN();
  if (!std::isfinite(static_cast<double>(grid_map.resolution)) ||
      grid_map.resolution <= 0.0f || !std::isfinite(grid_map.origin_x) ||
      !std::isfinite(grid_map.origin_y) || !std::isfinite(point.first) ||
      !std::isfinite(point.second))
  {
    return false;
  }

  const double converted_x = grid_map.origin_x +
    point.first * static_cast<double>(grid_map.resolution);
  const double converted_y = grid_map.origin_y +
    point.second * static_cast<double>(grid_map.resolution);
  if (!std::isfinite(converted_x) || !std::isfinite(converted_y))
    return false;

  wx = converted_x;
  wy = converted_y;
  return true;
}

// Keep interpolation in continuous grid coordinates.  In particular, do not
// round samples back to cells: doing so would move a sub-cell endpoint and
// could make the published path differ from the path that Core validated.
bool countInterpolatedPathPoints(
  const PathSolution& solution, size_t max_points, size_t& point_count,
  std::string& error)
{
  point_count = 0;
  error.clear();
  if (max_points < 2)
  {
    error = "max_path_points must be at least 2";
    return false;
  }
  if (solution.turning_points_.size() > max_points - 2)
  {
    error = "path has more essential waypoints than max_path_points=" +
      std::to_string(max_points);
    return false;
  }
  if (solution.turning_points_.size() >
      std::numeric_limits<size_t>::max() - 2)
  {
    error = "path has too many essential waypoints to represent";
    return false;
  }

  std::vector<Point2d> projected;
  projected.reserve(solution.turning_points_.size() + 2);
  projected.emplace_back(solution.start_);
  for (const auto& point : solution.turning_points_)
  {
    projected.emplace_back(
      static_cast<double>(point.first), static_cast<double>(point.second));
  }
  projected.emplace_back(solution.goal_);

  // Each segment contributes its first point and the final goal is appended
  // once.  Check the complete count before allocating the interpolation
  // vector; no rejected path can trigger an oversized reserve().
  size_t count = 1;
  for (size_t i = 0; i + 1 < projected.size(); ++i)
  {
    const auto& previous = projected[i];
    const auto& next = projected[i + 1];
    const double distance = std::hypot(
      previous.first - next.first, previous.second - next.second);
    if (!std::isfinite(distance))
    {
      error = "path contains a non-finite segment length";
      return false;
    }
    const double rounded_count = std::max(1.0, std::ceil(distance));
    if (rounded_count > static_cast<double>(max_points) ||
        rounded_count >= static_cast<double>(std::numeric_limits<size_t>::max()))
    {
      error = "path interpolation exceeds max_path_points=" +
        std::to_string(max_points);
      return false;
    }
    const size_t segment_count = static_cast<size_t>(rounded_count);
    if (segment_count > max_points - count)
    {
      error = "path interpolation exceeds max_path_points=" +
        std::to_string(max_points);
      return false;
    }
    count += segment_count;
  }

  point_count = count;
  return true;
}

bool interpolateProjectedPath(
  const PathSolution& solution, size_t max_points,
  std::vector<Point2d>& interpolated, std::string& error)
{
  interpolated.clear();
  size_t point_count = 0;
  if (!countInterpolatedPathPoints(
      solution, max_points, point_count, error))
  {
    return false;
  }

  std::vector<Point2d> projected;
  projected.reserve(solution.turning_points_.size() + 2);
  projected.emplace_back(solution.start_);
  for (const auto& point : solution.turning_points_)
  {
    projected.emplace_back(
      static_cast<double>(point.first), static_cast<double>(point.second));
  }
  projected.emplace_back(solution.goal_);

  interpolated.reserve(point_count);
  for (size_t i = 0; i + 1 < projected.size(); ++i)
  {
    const auto& previous = projected[i];
    const auto& next = projected[i + 1];
    const double distance = std::hypot(
      previous.first - next.first, previous.second - next.second);
    // Map dimensions are bounded by signed int in Core, so a finite segment
    // has a representable sample count on supported platforms.  Clamp the
    // degenerate case to one sample to avoid a zero divisor.
    const size_t count = static_cast<size_t>(std::max(
      1.0, std::ceil(distance)));
    for (size_t j = 0; j < count; ++j)
    {
      const double fraction = static_cast<double>(j) /
        static_cast<double>(count);
      interpolated.emplace_back(
        previous.first + (next.first - previous.first) * fraction,
        previous.second + (next.second - previous.second) * fraction);
    }
  }
  // Appending the final point explicitly preserves the continuous goal as the
  // last path pose instead of relying on a rounded interpolation sample.
  interpolated.emplace_back(projected.back());
  return true;
}

bool checkedAdd(size_t& total, size_t addition, size_t limit)
{
  if (addition > limit || total > limit - addition)
    return false;
  total += addition;
  return true;
}

bool estimatePathResponseBytes(
  size_t point_count, size_t frame_id_bytes, size_t& bytes)
{
  if (frame_id_bytes > std::numeric_limits<size_t>::max() - kPoseBaseBytes ||
      frame_id_bytes > std::numeric_limits<size_t>::max() - kPathBaseBytes)
    return false;
  const size_t bytes_per_pose = kPoseBaseBytes + frame_id_bytes;
  if (point_count >
      (std::numeric_limits<size_t>::max() - kPathBaseBytes - frame_id_bytes) /
        bytes_per_pose)
  {
    return false;
  }
  bytes = kPathBaseBytes + frame_id_bytes + point_count * bytes_per_pose;
  return true;
}

bool appendResponseNotice(
  std::string& message, const std::string& notice, size_t& response_bytes,
  size_t max_response_bytes)
{
  // Keep diagnostics bounded independently of caller-controlled frame IDs.
  if (message.size() >= kMaxDiagnosticBytes)
    return false;
  if (response_bytes >= max_response_bytes)
    return false;

  const size_t diagnostic_room = kMaxDiagnosticBytes - message.size();
  const size_t response_room = max_response_bytes - response_bytes;
  size_t room = std::min(diagnostic_room, response_room);
  const size_t separator_bytes = message.empty() ? 0u : std::min<size_t>(2, room);
  if (!message.empty() && room <= separator_bytes)
    return false;
  room -= separator_bytes;
  const size_t notice_bytes = std::min(room, notice.size());
  if (notice_bytes == 0)
    return false;

  if (separator_bytes != 0)
    message.append("; ", separator_bytes);
  message.append(notice, 0, notice_bytes);
  response_bytes += separator_bytes + notice_bytes;
  return notice_bytes == notice.size();
}

void setBoundedExceptionMessage(std::string& message, const char* detail)
{
  constexpr const char* prefix = "Raystar request failed with exception: ";
  message.assign(prefix);
  if (message.size() >= kMaxDiagnosticBytes || detail == nullptr)
  {
    message.resize(std::min(message.size(), kMaxDiagnosticBytes));
    return;
  }

  size_t detail_bytes = 0;
  const size_t remaining = kMaxDiagnosticBytes - message.size();
  while (detail_bytes < remaining && detail[detail_bytes] != '\0')
    ++detail_bytes;
  message.append(detail, detail_bytes);
}

uint32_t boundedUint32(size_t value)
{
  return static_cast<uint32_t>(std::min(
    value, static_cast<size_t>(std::numeric_limits<uint32_t>::max())));
}

uint64_t boundedUint64(size_t value)
{
  if constexpr (sizeof(size_t) > sizeof(uint64_t))
  {
    return static_cast<uint64_t>(std::min(
      value, static_cast<size_t>(std::numeric_limits<uint64_t>::max())));
  }
  return static_cast<uint64_t>(value);
}

uint16_t planningLimitMask(PlanningLimitReached limit)
{
  switch (limit)
  {
    case PlanningLimitReached::none:
      return PlanningResultInfo::LIMIT_NONE;
    case PlanningLimitReached::timeout:
      return PlanningResultInfo::LIMIT_TIMEOUT;
    case PlanningLimitReached::max_nodes:
      return PlanningResultInfo::LIMIT_MAX_NODES;
    case PlanningLimitReached::max_path_points:
      return PlanningResultInfo::LIMIT_MAX_PATH_POINTS;
    case PlanningLimitReached::cancelled:
      return PlanningResultInfo::LIMIT_CANCELLED;
  }
  return PlanningResultInfo::LIMIT_NONE;
}

void markPathOutputIncomplete(
  PlanningResultInfo& info, uint16_t limit = PlanningResultInfo::LIMIT_NONE)
{
  info.output_complete = false;
  info.limits_reached = static_cast<uint16_t>(info.limits_reached | limit);
  // Search interruption and terminal/error outcomes take precedence in the
  // high-level status. The independent output flag and limit bit still retain
  // an additional ROS serialization truncation.
  if (info.status == PlanningResultInfo::STATUS_COMPLETE ||
      info.status == PlanningResultInfo::STATUS_FEWER_PATHS ||
      info.status == PlanningResultInfo::STATUS_NO_PATH)
  {
    info.status = PlanningResultInfo::STATUS_PARTIAL_OUTPUT;
  }
}

template<typename ResponseT>
void resetPlanningResponse(ResponseT& response)
{
  response.success = false;
  response.result_info = PlanningResultInfo{};
  response.message.clear();
  response.path_results.clear();
  response.debug_nodes.clear();
}

template<typename ResponseT>
void setRequestedPathCount(ResponseT& response, int32_t requested)
{
  response.result_info.requested_path_count = requested > 0 ?
    static_cast<uint32_t>(requested) : 0u;
}

template<typename ResponseT>
void initializePlanningResponse(
  ResponseT& response, int32_t requested,
  const raystar_interfaces::MapId& map_id, bool debug_requested)
{
  resetPlanningResponse(response);
  response.result_info.map_id = map_id;
  response.result_info.debug_requested = debug_requested;
  response.result_info.debug_output_complete = !debug_requested;
  setRequestedPathCount(response, requested);
}

template<typename ResponseT>
void markCancelled(ResponseT& response)
{
  response.success = false;
  response.path_results.clear();
  response.debug_nodes.clear();
  response.result_info.status = PlanningResultInfo::STATUS_CANCELLED;
  response.result_info.limits_reached = static_cast<uint16_t>(
    response.result_info.limits_reached |
    PlanningResultInfo::LIMIT_CANCELLED);
  response.result_info.request_satisfied = false;
  response.result_info.search_complete = false;
  response.result_info.output_complete = false;
  response.result_info.debug_output_complete = false;
  response.result_info.returned_path_count = 0;
}

template<typename ResponseT>
void markFailed(ResponseT& response)
{
  response.success = false;
  response.path_results.clear();
  response.debug_nodes.clear();
  response.result_info.status = PlanningResultInfo::STATUS_FAILED;
  response.result_info.request_satisfied = false;
  response.result_info.search_complete = false;
  response.result_info.output_complete = false;
  response.result_info.debug_output_complete = false;
  response.result_info.returned_path_count = 0;
}

class PlanningSlotRelease
{
public:
  explicit PlanningSlotRelease(std::atomic<bool>& busy) : busy_(busy) {}
  PlanningSlotRelease(const PlanningSlotRelease&) = delete;
  PlanningSlotRelease& operator=(const PlanningSlotRelease&) = delete;
  ~PlanningSlotRelease()
  {
    busy_.store(false, std::memory_order_release);
  }

private:
  std::atomic<bool>& busy_;
};

class SearchStateRelease
{
public:
  explicit SearchStateRelease(RaystarCore& core) : core_(core) {}
  SearchStateRelease(const SearchStateRelease&) = delete;
  SearchStateRelease& operator=(const SearchStateRelease&) = delete;
  ~SearchStateRelease()
  {
    core_.resetSearchState();
  }

private:
  RaystarCore& core_;
};

template<typename Callback>
class ScopeExit
{
public:
  explicit ScopeExit(Callback callback)
    : callback_(std::move(callback)) {}
  ScopeExit(const ScopeExit&) = delete;
  ScopeExit& operator=(const ScopeExit&) = delete;
  ~ScopeExit() noexcept
  {
    try
    {
      callback_();
    }
    catch (...)
    {
    }
  }

private:
  Callback callback_;
};

}  // namespace

RaystarNode::RaystarNode(const rclcpp::NodeOptions& options)
  : Node("raystar", options)
{
  // Register before declaration so command-line/YAML overrides are checked by
  // the same validator as runtime updates.  The callback has no side effects:
  // rclcpp owns the all-or-none parameter commit after every callback accepts.
  parameter_validation_callback_ = add_on_set_parameters_callback(
    validateParameterChanges);
  for (const auto& spec : kIntegerParameterSpecs)
  {
    declare_parameter<int64_t>(
      spec.name, spec.default_value, makeParameterDescriptor(spec));
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
  declare_parameter<bool>(
    "enable_legacy_map_service", true, legacy_service_descriptor);

  RequestConfiguration configuration;
  std::string configuration_error;
  if (!loadRequestConfiguration(configuration, configuration_error))
  {
    throw rclcpp::exceptions::InvalidParameterValueException(
      "Invalid Raystar server configuration: " + configuration_error);
  }

  if (get_parameter("enable_legacy_map_service").as_bool())
  {
    service_ = create_service<raystar_interfaces::srv::GetRaystarPaths>(
      "~/get_raystar_paths",
      std::bind(&RaystarNode::handleService, this,
        std::placeholders::_1, std::placeholders::_2));
    RCLCPP_WARN(get_logger(),
      "Legacy full-map Service is enabled; prefer the cached-map Action and "
      "disable enable_legacy_map_service in resource-constrained deployments");
  }

  map_status_publisher_ =
    create_publisher<raystar_interfaces::msg::MapStatus>(
      "~/map_status", rclcpp::QoS(rclcpp::KeepLast(1)).transient_local());
  const std::string map_topic = get_parameter("map_topic").as_string();
  map_subscription_ = create_subscription<nav_msgs::msg::OccupancyGrid>(
    map_topic, rclcpp::QoS(rclcpp::KeepLast(1)).transient_local().reliable(),
    [this](nav_msgs::msg::OccupancyGrid::ConstSharedPtr map) {
      handleMap(std::move(map));
    });

  action_server_ = rclcpp_action::create_server<PlanAction>(
    this, "~/plan_paths",
    std::bind(&RaystarNode::handleActionGoal, this,
      std::placeholders::_1, std::placeholders::_2),
    std::bind(&RaystarNode::handleActionCancel, this,
      std::placeholders::_1),
    std::bind(&RaystarNode::handleActionAccepted, this,
      std::placeholders::_1));

  auto latched = rclcpp::QoS(rclcpp::KeepLast(1)).transient_local();
  non_homotopic_pub_ = create_publisher<visualization_msgs::msg::MarkerArray>(
    "~/non_homotopic_paths", latched);
  poly_obstacle_pub_ = create_publisher<visualization_msgs::msg::MarkerArray>(
    "~/poly_obstacles", latched);
  debug_tree_pub_ = create_publisher<visualization_msgs::msg::MarkerArray>(
    "~/debug_tree", latched);
  cdt_pub_ = create_publisher<visualization_msgs::msg::MarkerArray>(
    "~/cdt", latched);

  const int64_t path_republish_period_ms =
    get_parameter("path_visualization_republish_period_ms").as_int();
  if (path_republish_period_ms > 0)
  {
    path_visualization_timer_ = create_wall_timer(
      std::chrono::milliseconds(path_republish_period_ms),
      std::bind(&RaystarNode::republishCachedPathVisualization, this));
  }

  RCLCPP_INFO(get_logger(),
    "Raystar node initialized: occupied_threshold=%d "
    "max_k=%d max_nodes=%zu "
    "planning_timeout_ms=%lld max_map_cells=%zu max_map_bytes=%zu "
    "max_path_points=%zu max_debug_nodes=%zu max_response_bytes=%zu "
    "path_visualization_republish_period_ms=%lld map_topic=%s",
    configuration.occupied_threshold,
    configuration.planning.max_k, configuration.planning.max_nodes,
    static_cast<long long>(configuration.planning.planning_timeout.count()),
    configuration.planning.max_map_cells,
    configuration.planning.max_map_bytes,
    configuration.planning.max_path_points,
    configuration.planning.max_debug_nodes,
    configuration.planning.max_response_bytes,
    static_cast<long long>(path_republish_period_ms), map_topic.c_str());

  // Start the persistent capacity-one Action worker only after every other
  // constructor operation that can throw has completed.  Accepted/cancel
  // callbacks only update bounded state and notify this worker.
  action_worker_ = std::thread(&RaystarNode::actionWorkerLoop, this);
}

RaystarNode::~RaystarNode()
{
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

bool RaystarNode::loadRequestConfiguration(
  RequestConfiguration& configuration, std::string& error) const
{
  const auto parameters = get_parameters(
    {"max_k", "max_nodes", "planning_timeout_ms", "max_map_cells",
     "max_map_bytes", "max_path_points", "max_debug_nodes",
     "max_response_bytes", "occupied_threshold"});
  if (parameters.size() != 9)
  {
    error = "could not read the complete ROS server parameter set";
    return false;
  }
  const auto validation = validateParameterChanges(parameters);
  if (!validation.successful)
  {
    error = validation.reason;
    return false;
  }

  const int64_t max_k = parameters[0].as_int();
  const int64_t max_nodes = parameters[1].as_int();
  const int64_t planning_timeout_ms = parameters[2].as_int();
  const int64_t max_map_cells = parameters[3].as_int();
  const int64_t max_map_bytes = parameters[4].as_int();
  const int64_t max_path_points = parameters[5].as_int();
  const int64_t max_debug_nodes = parameters[6].as_int();
  const int64_t max_response_bytes = parameters[7].as_int();
  const int64_t occupied_threshold = parameters[8].as_int();
  configuration = RequestConfiguration{};
  configuration.occupied_threshold = static_cast<int>(occupied_threshold);
  configuration.planning.max_k = static_cast<int>(max_k);
  configuration.planning.max_nodes = static_cast<size_t>(max_nodes);
  configuration.planning.planning_timeout =
    std::chrono::milliseconds(planning_timeout_ms);
  configuration.planning.max_map_cells = static_cast<size_t>(max_map_cells);
  configuration.planning.max_map_bytes = static_cast<size_t>(max_map_bytes);
  configuration.planning.max_path_points = static_cast<size_t>(max_path_points);
  configuration.planning.max_debug_nodes = static_cast<size_t>(max_debug_nodes);
  configuration.planning.max_response_bytes =
    static_cast<size_t>(max_response_bytes);
  error.clear();
  return true;
}

void RaystarNode::handleMap(
  nav_msgs::msg::OccupancyGrid::ConstSharedPtr map)
{
  if (!map || shutting_down_.load(std::memory_order_acquire))
    return;

  try
  {
    RequestConfiguration configuration;
    std::string error;
    if (!loadRequestConfiguration(configuration, error))
    {
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
        *map, false, configuration, stop_requested, validation_map, error))
    {
      RCLCPP_WARN(get_logger(), "Rejected map cache update: %s", error.c_str());
      return;
    }

    const raystar_interfaces::MapId id =
      raystar_interfaces::computeMapId(*map);
    std::uint64_t generation = 0;
    {
      std::lock_guard<std::mutex> lock(map_cache_mutex_);
      generation = cached_map_.generation ==
        std::numeric_limits<std::uint64_t>::max() ?
        1U : cached_map_.generation + 1U;
      cached_map_.message = map;
      cached_map_.id = id;
      cached_map_.generation = generation;
    }
    raystar_interfaces::msg::MapStatus status;
    status.header = map->header;
    status.map_id = id;
    status.generation = generation;
    status.width = map->info.width;
    status.height = map->info.height;
    status.resolution = map->info.resolution;
    map_status_publisher_->publish(status);
    RCLCPP_INFO(get_logger(),
      "Cached admitted OccupancyGrid generation=%llu",
      static_cast<unsigned long long>(generation));
  }
  catch (const std::exception& exception)
  {
    RCLCPP_ERROR(get_logger(), "Could not cache OccupancyGrid: %s",
      exception.what());
  }
  catch (...)
  {
    RCLCPP_ERROR(get_logger(),
      "Could not cache OccupancyGrid because of an unknown exception");
  }
}

bool RaystarNode::resolveCachedMap(
  const raystar_interfaces::MapId& requested_id,
  nav_msgs::msg::OccupancyGrid::ConstSharedPtr& map,
  std::string& error) const
{
  map.reset();
  if (raystar_interfaces::isZeroMapId(requested_id))
  {
    error = "Invalid map_id: cached-map Action requires a non-zero map identity";
    return false;
  }

  std::lock_guard<std::mutex> lock(map_cache_mutex_);
  if (!cached_map_.message)
  {
    error = "No admitted OccupancyGrid is cached; publish the configured map_topic first";
    return false;
  }
  if (!raystar_interfaces::mapIdsEqual(requested_id, cached_map_.id))
  {
    error =
      "Requested map_id is not the current cached OccupancyGrid; refresh the map snapshot";
    return false;
  }
  map = cached_map_.message;
  error.clear();
  return true;
}

bool RaystarNode::occupancyGridToBinaryMap(
  const nav_msgs::msg::OccupancyGrid& grid, bool allow_unknown,
  const RequestConfiguration& configuration,
  const StopPredicate& stop_requested,
  GridMap& output, std::string& error) const
{
  output = GridMap{};
  error.clear();
  const size_t width = static_cast<size_t>(grid.info.width);
  const size_t height = static_cast<size_t>(grid.info.height);

  if (grid.header.frame_id.empty())
  {
    error = "Invalid map: header.frame_id must not be empty";
    return false;
  }

  if (!std::isfinite(static_cast<double>(grid.info.resolution)) ||
      grid.info.resolution <= 0.0f)
  {
    error = "Invalid map: resolution must be finite and greater than zero";
    return false;
  }

  const auto& origin = grid.info.origin;
  if (!std::isfinite(origin.position.x) ||
      !std::isfinite(origin.position.y) ||
      !std::isfinite(origin.position.z))
  {
    error = "Invalid map: origin position must contain only finite coordinates";
    return false;
  }
  if (origin.position.z != 0.0)
  {
    error = "Invalid map: origin z coordinate must be zero for 2D planning";
    return false;
  }

  const auto& orientation = origin.orientation;
  if (!std::isfinite(orientation.x) || !std::isfinite(orientation.y) ||
      !std::isfinite(orientation.z) || !std::isfinite(orientation.w))
  {
    error = "Invalid map: origin orientation must contain only finite values";
    return false;
  }
  const double quaternion_norm_squared =
    orientation.x * orientation.x + orientation.y * orientation.y +
    orientation.z * orientation.z + orientation.w * orientation.w;
  if (!std::isfinite(quaternion_norm_squared) ||
      std::abs(quaternion_norm_squared - 1.0) > kQuaternionNormTolerance)
  {
    error = "Invalid map: origin orientation must be a unit quaternion";
    return false;
  }
  if (orientation.x != 0.0 || orientation.y != 0.0 || orientation.z != 0.0)
  {
    error = "Invalid map: origin orientation must be the identity rotation; "
      "rotated or tilted occupancy grids are not supported";
    return false;
  }

  // Reuse Core's allocation-free admission check so ROS and direct callers
  // cannot drift on overflow, signed-indexing, or byte-budget semantics.
  MapResourceEstimate map_estimate;
  if (!validateMapResourceBudget(
      width, height, grid.data.size(), configuration.planning,
      map_estimate, error))
  {
    return false;
  }
  const size_t cell_count = map_estimate.cell_count;

  const double max_world_x = origin.position.x +
    static_cast<double>(width) * static_cast<double>(grid.info.resolution);
  const double max_world_y = origin.position.y +
    static_cast<double>(height) * static_cast<double>(grid.info.resolution);
  if (!std::isfinite(max_world_x) || !std::isfinite(max_world_y))
  {
    error = "Invalid map: world-coordinate extent is not finite";
    return false;
  }

  GridMap map;
  map.width = grid.info.width;
  map.height = grid.info.height;
  map.resolution = grid.info.resolution;
  map.origin_x = grid.info.origin.position.x;
  map.origin_y = grid.info.origin.position.y;

  map.data.resize(cell_count);
  for (size_t i = 0; i < cell_count; ++i)
  {
    if ((i & 0xfffu) == 0u && stop_requested && stop_requested())
    {
      error = "Planning was cancelled while converting the occupancy grid";
      return false;
    }
    const int value = static_cast<int>(grid.data[i]);
    if (value < -1 || value > 100)
    {
      const size_t x = i % width;
      const size_t y = i / width;
      error = "Invalid map: occupancy value " + std::to_string(value) +
        " at cell (" + std::to_string(x) + ", " + std::to_string(y) +
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
  std::shared_ptr<raystar_interfaces::srv::GetRaystarPaths::Response> response)
{
  resetPlanningResponse(*response);
  setRequestedPathCount(*response, request ? request->k : 0);
  if (!request)
  {
    response->result_info.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
    response->message = "Raystar Service request is unavailable";
    return;
  }
  bool expected_idle = false;
  if (!planning_busy_.compare_exchange_strong(
      expected_idle, true, std::memory_order_acq_rel))
  {
    response->result_info.status = PlanningResultInfo::STATUS_BUSY;
    response->message =
      "Raystar planner is busy with another Service request or Action goal";
    return;
  }

  PlanningSlotRelease release_slot(planning_busy_);
  const auto map_id = raystar_interfaces::computeMapId(request->map);
  executePlanning(*request, *response, request->map, map_id, [this]() {
    return shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok();
  });
}

rclcpp_action::GoalResponse RaystarNode::handleActionGoal(
  const rclcpp_action::GoalUUID& uuid,
  std::shared_ptr<const PlanAction::Goal> goal)
{
  if (!goal || shutting_down_.load(std::memory_order_acquire) ||
      !rclcpp::ok())
  {
    return rclcpp_action::GoalResponse::REJECT;
  }

  std::shared_ptr<std::atomic<bool>> cancel_requested;
  try
  {
    cancel_requested = std::make_shared<std::atomic<bool>>(false);
  }
  catch (...)
  {
    RCLCPP_ERROR(get_logger(),
      "Rejecting Action goal: could not allocate cancellation state");
    return rclcpp_action::GoalResponse::REJECT;
  }

  bool expected_idle = false;
  if (!planning_busy_.compare_exchange_strong(
      expected_idle, true, std::memory_order_acq_rel))
  {
    RCLCPP_WARN(get_logger(),
      "Rejecting Action goal because the capacity-one planner is busy");
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
  const std::shared_ptr<PlanGoalHandle> goal_handle)
{
  if (!goal_handle)
    return rclcpp_action::CancelResponse::REJECT;

  std::lock_guard<std::mutex> lock(action_state_mutex_);
  if (!active_goal_reserved_ ||
      active_goal_id_ != goal_handle->get_goal_id() ||
      !active_goal_cancel_)
  {
    return rclcpp_action::CancelResponse::REJECT;
  }
  // Core polls this per-goal flag directly.  rclcpp_action transitions the
  // goal to CANCELING immediately after this callback returns ACCEPT; the
  // worker waits for that transition before publishing its canceled result.
  active_goal_cancel_->store(true, std::memory_order_release);
  return rclcpp_action::CancelResponse::ACCEPT;
}

void RaystarNode::handleActionAccepted(
  const std::shared_ptr<PlanGoalHandle> goal_handle)
{
  std::shared_ptr<std::atomic<bool>> cancel_requested;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (goal_handle && active_goal_reserved_ &&
        active_goal_id_ == goal_handle->get_goal_id())
    {
      cancel_requested = active_goal_cancel_;
    }
  }

  if (!goal_handle || !cancel_requested)
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    active_goal_reserved_ = false;
    active_goal_cancel_.reset();
    planning_busy_.store(false, std::memory_order_release);
    return;
  }

  bool queued = false;
  {
    std::lock_guard<std::mutex> lock(action_worker_mutex_);
    if (!stop_action_worker_ && !pending_action_job_)
    {
      pending_action_job_.emplace(
        ActionJob{goal_handle, cancel_requested});
      queued = true;
    }
  }
  if (queued)
  {
    action_worker_cv_.notify_one();
    return;
  }

  // Release admission before doing any best-effort result allocation. This
  // branch is only reachable during shutdown or an internal invariant
  // violation; even an allocation failure must not leave the planner busy.
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (active_goal_reserved_ &&
        active_goal_id_ == goal_handle->get_goal_id())
    {
      active_goal_reserved_ = false;
      active_goal_cancel_.reset();
    }
    planning_busy_.store(false, std::memory_order_release);
  }
  try
  {
    auto result = std::make_shared<PlanAction::Result>();
    const auto rejected_goal =
      goal_handle ? goal_handle->get_goal() : nullptr;
    if (rejected_goal)
    {
      initializePlanningResponse(
        *result, rejected_goal->k, rejected_goal->map_id,
        rejected_goal->include_debug);
    }
    else
    {
      resetPlanningResponse(*result);
    }
    result->result_info.status = PlanningResultInfo::STATUS_FAILED;
    result->message = "Raystar Action worker is unavailable";
    goal_handle->abort(result);
  }
  catch (...)
  {
  }
}

void RaystarNode::actionWorkerLoop() noexcept
{
  while (true)
  {
    ActionJob job;
    {
      std::unique_lock<std::mutex> lock(action_worker_mutex_);
      action_worker_cv_.wait(lock, [this]() {
        return stop_action_worker_ || pending_action_job_.has_value();
      });
      if (stop_action_worker_ && !pending_action_job_)
        return;
      job = std::move(*pending_action_job_);
      pending_action_job_.reset();
    }
    executeAction(job.goal_handle, job.cancel_requested);
  }
}

void RaystarNode::executeAction(
  const std::shared_ptr<PlanGoalHandle> goal_handle,
  const std::shared_ptr<std::atomic<bool>>& cancel_requested) noexcept
{
  enum class TerminalState
  {
    succeeded,
    aborted,
    canceled
  };

  std::shared_ptr<PlanAction::Result> result;
  TerminalState terminal_state = TerminalState::aborted;
  const auto goal_is_canceling = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try
    {
      return goal_handle && goal_handle->is_canceling();
    }
    catch (...)
    {
      return false;
    }
  };
  const auto goal_is_active = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try
    {
      return goal_handle && goal_handle->is_active();
    }
    catch (...)
    {
      return false;
    }
  };

  try
  {
    result = std::make_shared<PlanAction::Result>();
    resetPlanningResponse(*result);
    const std::weak_ptr<PlanGoalHandle> weak_goal_handle(goal_handle);
    const StopPredicate stop_requested =
      [this, weak_goal_handle, cancel_requested]() {
        if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
          return true;
        const auto handle = weak_goal_handle.lock();
        if (!handle)
          return true;
        // Stop Core as soon as our cancel callback accepts the request.  The
        // worker separately waits for rclcpp_action's subsequent CANCELING
        // transition before it publishes the terminal result.
        return cancel_requested->load(std::memory_order_acquire);
      };

    const auto goal = goal_handle->get_goal();
    if (goal)
    {
      nav_msgs::msg::OccupancyGrid::ConstSharedPtr cached_map;
      std::string map_error;
      if (!resolveCachedMap(goal->map_id, cached_map, map_error))
      {
        initializePlanningResponse(
          *result, goal->k, goal->map_id, goal->include_debug);
        result->result_info.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
        result->message = std::move(map_error);
      }
      else
      {
        executePlanning(
          *goal, *result, *cached_map, goal->map_id, stop_requested);
      }
    }
    else
    {
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->message = "Raystar Action goal data is unavailable";
    }

    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
    {
      markFailed(*result);
      if (result->message.empty())
        result->message = "Planning stopped because the Raystar node is shutting down";
      terminal_state = TerminalState::aborted;
    }
    else if (goal_is_canceling())
    {
      markCancelled(*result);
      if (result->message.empty())
        result->message = "Planning was cancelled";
      terminal_state = TerminalState::canceled;
    }
    else if (result->success)
    {
      terminal_state = TerminalState::succeeded;
    }
    else
    {
      terminal_state = TerminalState::aborted;
    }
  }
  catch (const std::exception& exception)
  {
    try
    {
      const uint32_t requested_path_count = result ?
        result->result_info.requested_path_count : 0u;
      if (!result)
        result = std::make_shared<PlanAction::Result>();
      resetPlanningResponse(*result);
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->result_info.requested_path_count = requested_path_count;
      setBoundedExceptionMessage(result->message, exception.what());
    }
    catch (...)
    {
    }
    terminal_state = goal_is_canceling() ?
      TerminalState::canceled : TerminalState::aborted;
    RCLCPP_ERROR(get_logger(), "Raystar Action worker failed: %s",
      exception.what());
  }
  catch (...)
  {
    try
    {
      const uint32_t requested_path_count = result ?
        result->result_info.requested_path_count : 0u;
      if (!result)
        result = std::make_shared<PlanAction::Result>();
      resetPlanningResponse(*result);
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->result_info.requested_path_count = requested_path_count;
      result->message = "Raystar Action worker failed with an unknown exception";
    }
    catch (...)
    {
    }
    terminal_state = goal_is_canceling() ?
      TerminalState::canceled : TerminalState::aborted;
    RCLCPP_ERROR(get_logger(),
      "Raystar Action worker failed with an unknown exception");
  }

  // Linearize completion against handleActionCancel().  If the cancel
  // callback stored its per-goal flag first, retain the reservation while
  // rclcpp_action performs its internal transition after that callback
  // returns.  Otherwise clear the reservation under the same mutex, causing a
  // later cancel request to be rejected before this normal terminal result.
  bool explicit_cancel = false;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    explicit_cancel =
      cancel_requested->load(std::memory_order_acquire) &&
      !shutting_down_.load(std::memory_order_acquire);
    if (!explicit_cancel)
    {
      if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
        terminal_state = TerminalState::aborted;
      if (active_goal_reserved_ &&
          active_goal_id_ == goal_handle->get_goal_id())
      {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
      planning_busy_.store(false, std::memory_order_release);
    }
  }

  // An accepted cancel must retain the capacity-one reservation through both
  // Marker/cache cleanup and the Action terminal transition.  Releasing busy
  // earlier would allow a new request to publish a fresh snapshot which this
  // older goal could then erase in its normal or abort-to-cancel fallback.
  ScopeExit explicit_cancel_release([&]() noexcept {
    if (!explicit_cancel)
      return;
    try
    {
      std::lock_guard<std::mutex> lock(action_state_mutex_);
      if (active_goal_reserved_ &&
          active_goal_id_ == goal_handle->get_goal_id())
      {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
    }
    catch (...)
    {
      // Admission still has to be released if middleware state is already
      // tearing down; the node destructor owns the remaining cleanup.
    }
    planning_busy_.store(false, std::memory_order_release);
  });

  if (explicit_cancel)
  {
    terminal_state = TerminalState::canceled;
    const auto transition_deadline =
      std::chrono::steady_clock::now() + std::chrono::seconds(1);
    while (rclcpp::ok() && goal_is_active() && !goal_is_canceling() &&
      std::chrono::steady_clock::now() < transition_deadline)
    {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (!goal_is_canceling())
    {
      terminal_state = TerminalState::aborted;
      if (result)
      {
        const uint32_t requested_path_count =
          result->result_info.requested_path_count;
        resetPlanningResponse(*result);
        result->result_info.requested_path_count = requested_path_count;
        markFailed(*result);
        result->message =
          "Cancellation was accepted but the Action state transition did not complete";
      }
    }

  }

  // Allocation failure is the only case in which no typed Action result can
  // be published.  The reservation has still been released and destruction
  // of the goal handle will force a terminal state inside rclcpp_action.
  if (!result)
    return;

  if (terminal_state == TerminalState::canceled)
  {
    markCancelled(*result);
    if (result->message.empty())
      result->message = "Planning was cancelled";
    std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
    clearVisualizationsLocked();
  }

  try
  {
    if (terminal_state == TerminalState::canceled)
      goal_handle->canceled(result);
    else if (terminal_state == TerminalState::succeeded)
      goal_handle->succeed(result);
    else
      goal_handle->abort(result);
  }
  catch (const std::exception& exception)
  {
    // The bounded cancel-transition fallback can race with the transition at
    // its deadline.  If abort() lost that race, finish with the now-valid
    // canceled transition instead of leaving the goal non-terminal.
    if (terminal_state != TerminalState::canceled &&
        goal_is_canceling())
    {
      try
      {
        markCancelled(*result);
        result->message = "Planning was cancelled";
        {
          std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
          clearVisualizationsLocked();
        }
        goal_handle->canceled(result);
        return;
      }
      catch (...)
      {
      }
    }
    RCLCPP_ERROR(get_logger(),
      "Could not publish Raystar Action terminal result: %s",
      exception.what());
  }
  catch (...)
  {
    RCLCPP_ERROR(get_logger(),
      "Could not publish Raystar Action terminal result");
  }
}

template<typename RequestT, typename ResponseT>
void RaystarNode::executePlanning(
  const RequestT& request_value, ResponseT& response_value,
  const nav_msgs::msg::OccupancyGrid& grid,
  const raystar_interfaces::MapId& map_id,
  const StopPredicate& stop_requested) noexcept try
{
  // The Core search tree and retained visualization cache are one stateful
  // unit.  The timer uses try_lock(), so it skips a tick instead of delaying
  // planning or Action cancellation while this lock is held.
  std::unique_lock<std::mutex> planner_lock(planner_cache_mutex_);
  const auto* request = &request_value;
  auto* response = &response_value;
  initializePlanningResponse(
    *response, request->k, map_id, request->include_debug);
  if (stop_requested && stop_requested())
  {
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
  if (!loadRequestConfiguration(configuration, configuration_error))
  {
    response->success = false;
    response->result_info.status =
      PlanningResultInfo::STATUS_INVALID_CONFIGURATION;
    response->message =
      "Invalid Raystar server configuration: " + configuration_error;
    return;
  }
  configuration.planning.cancel_requested = stop_requested;
  if (grid.header.frame_id.size() > kMaxFrameIdBytes ||
      request->start.header.frame_id.size() > kMaxFrameIdBytes ||
      request->goal.header.frame_id.size() > kMaxFrameIdBytes)
  {
    response->success = false;
    response->message =
      "Invalid request: frame_id must be at most 256 bytes";
    return;
  }
  if (grid.header.frame_id.empty())
  {
    response->success = false;
    response->message = "Invalid map: header.frame_id must not be empty";
    return;
  }
  if (request->k <= 0)
  {
    response->success = false;
    response->message = "Invalid K: must be greater than zero";
    return;
  }
  if (request->k > configuration.planning.max_k)
  {
    response->success = false;
    response->message = "Invalid K: requested " + std::to_string(request->k) +
      " exceeds max_k=" + std::to_string(configuration.planning.max_k);
    return;
  }

  // Validate the small, scalar endpoint fields before copying the map.  This
  // keeps an otherwise rejected request from paying the binary-map allocation
  // cost (the frame-length checks above provide the same guard for strings).
  std::string pose_error;
  if (!validatePlanarPose(
      request->start, "Start", grid.header.frame_id, pose_error))
  {
    response->success = false;
    response->message = pose_error;
    return;
  }
  if (!validatePlanarPose(
      request->goal, "Goal", grid.header.frame_id, pose_error))
  {
    response->success = false;
    response->message = pose_error;
    return;
  }

  GridMap work_map;
  std::string map_error;
  if (!occupancyGridToBinaryMap(
      grid, request->allow_unknown, configuration, stop_requested,
      work_map, map_error))
  {
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
      work_map, request->start.pose.position.x, request->start.pose.position.y,
      start_point)) {
    response->success = false;
    response->message = "Start position is outside the map";
    return;
  }
  if (!worldToContinuousMap(
      work_map, request->goal.pose.position.x, request->goal.pose.position.y,
      goal_point)) {
    response->success = false;
    response->message = "Goal position is outside the map";
    return;
  }

  const PlanEndpoint start_endpoint(
    {static_cast<int>(start_point.cell_x),
     static_cast<int>(start_point.cell_y)},
    {start_point.x, start_point.y});
  const PlanEndpoint goal_endpoint(
    {static_cast<int>(goal_point.cell_x),
     static_cast<int>(goal_point.cell_y)},
    {goal_point.x, goal_point.y});

  RCLCPP_INFO(get_logger(),
    "Planning: start=(%.6f,%.6f) cell=(%u,%u) goal=(%.6f,%.6f) "
    "cell=(%u,%u) K=%d allow_self_crossing=%s allow_unknown=%s "
    "occupied_threshold=%d",
    start_point.x, start_point.y, start_point.cell_x, start_point.cell_y,
    goal_point.x, goal_point.y, goal_point.cell_x, goal_point.cell_y,
    request->k,
    request->allow_self_crossing ? "true" : "false",
    request->allow_unknown ? "true" : "false",
    configuration.occupied_threshold);

  auto result = core_.plan(
    work_map, start_endpoint, goal_endpoint, request->k,
    request->allow_self_crossing, configuration.planning);
  SearchStateRelease search_state_release(core_);
  const auto& nodes = core_.getNodes();

  auto& result_info = response->result_info;
  result_info.found_path_count = boundedUint32(result.path_solutions.size());
  result_info.expanded_nodes = boundedUint64(result.expanded_nodes);
  result_info.map_time_ms = result.map_time_ms;
  result_info.plan_time_ms = result.plan_time_ms;
  result_info.limits_reached = planningLimitMask(result.limit_reached);
  switch (result.outcome)
  {
    case PlanningOutcome::complete:
      result_info.search_complete = true;
      result_info.status =
        result.path_solutions.size() >= static_cast<size_t>(request->k) ?
        PlanningResultInfo::STATUS_COMPLETE :
        PlanningResultInfo::STATUS_FEWER_PATHS;
      break;
    case PlanningOutcome::no_path:
      result_info.search_complete = true;
      result_info.status = PlanningResultInfo::STATUS_NO_PATH;
      break;
    case PlanningOutcome::limit_reached:
      result_info.search_complete = false;
      result_info.status =
        result.limit_reached == PlanningLimitReached::cancelled ?
        PlanningResultInfo::STATUS_CANCELLED :
        PlanningResultInfo::STATUS_PARTIAL_SEARCH;
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
    result.map_time_ms, result.plan_time_ms, result.path_solutions.size());

  if (stop_requested && stop_requested())
  {
    markCancelled(*response);
    response->message = result.message.empty() ?
      "Planning was cancelled" : result.message;
    return;
  }

  response->message = result.message;
  if (response->message.size() > kMaxDiagnosticBytes)
    response->message.resize(kMaxDiagnosticBytes);

  const double resolution = static_cast<double>(work_map.resolution);
  size_t response_bytes = kResponseBaseBytes;
  if (!checkedAdd(response_bytes, response->message.size(),
      configuration.planning.max_response_bytes) ||
      !checkedAdd(response_bytes, grid.header.frame_id.size(),
      configuration.planning.max_response_bytes))
  {
    response->success = false;
    markPathOutputIncomplete(
      result_info, PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
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
  std::vector<size_t> emitted_solution_indices;
  size_t omitted_path_count = 0;
  size_t emitted_path_points = 0;
  uint16_t path_output_limits = PlanningResultInfo::LIMIT_NONE;
  for (size_t solution_index = 0;
      solution_index < result.path_solutions.size(); ++solution_index)
  {
    if (stop_requested && stop_requested())
    {
      initializePlanningResponse(
        *response, request->k, map_id, request->include_debug);
      result_info.found_path_count = boundedUint32(result.path_solutions.size());
      result_info.expanded_nodes = boundedUint64(result.expanded_nodes);
      result_info.map_time_ms = result.map_time_ms;
      result_info.plan_time_ms = result.plan_time_ms;
      markCancelled(*response);
      response->message = "Planning was cancelled while building the response";
      return;
    }
    const auto& sol = result.path_solutions[solution_index];
    size_t point_count = 0;
    std::string path_error;
    if (!countInterpolatedPathPoints(
        sol, configuration.planning.max_path_points, point_count, path_error))
    {
      ++omitted_path_count;
      path_output_limits = static_cast<uint16_t>(
        path_output_limits | PlanningResultInfo::LIMIT_MAX_PATH_POINTS);
      appendResponseNotice(response->message,
        "A path was omitted: " + path_error, response_bytes,
        configuration.planning.max_response_bytes);
      continue;
    }
    if (point_count > configuration.planning.max_path_points -
        emitted_path_points)
    {
      omitted_path_count = result.path_solutions.size() -
        emitted_solution_indices.size();
      path_output_limits = static_cast<uint16_t>(
        path_output_limits | PlanningResultInfo::LIMIT_MAX_PATH_POINTS);
      appendResponseNotice(response->message,
        "Path output reached the per-response max_path_points=" +
        std::to_string(configuration.planning.max_path_points), response_bytes,
        configuration.planning.max_response_bytes);
      break;
    }

    size_t path_bytes = 0;
    if (!estimatePathResponseBytes(
        point_count, grid.header.frame_id.size(), path_bytes) ||
        !checkedAdd(response_bytes, path_bytes,
        configuration.planning.max_response_bytes))
    {
      ++omitted_path_count;
      path_output_limits = static_cast<uint16_t>(
        path_output_limits | PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
      appendResponseNotice(response->message,
        "Path output reached max_response_bytes=" +
        std::to_string(configuration.planning.max_response_bytes),
        response_bytes, configuration.planning.max_response_bytes);
      // Paths are ordered by Core's cost/topology ranking.  Once one no longer
      // fits, later paths cannot make the response smaller in a predictable
      // way, so stop before allocating another Path message.
      omitted_path_count = result.path_solutions.size() -
        emitted_solution_indices.size();
      break;
    }

    raystar_interfaces::msg::PathResult path_result;
    const size_t remaining_path_points =
      configuration.planning.max_path_points - emitted_path_points;
    if (!buildPathMsg(
        sol, work_map, grid.header.frame_id,
        remaining_path_points,
        path_result.path, path_error))
    {
      ++omitted_path_count;
      appendResponseNotice(response->message,
        "A path was omitted: " + path_error, response_bytes,
        configuration.planning.max_response_bytes);
      continue;
    }
    path_result.cost = sol.path_cost_ * resolution;
    response->path_results.emplace_back(std::move(path_result));
    emitted_solution_indices.push_back(solution_index);
    emitted_path_points += point_count;
  }

  if (omitted_path_count > 0)
  {
    appendResponseNotice(response->message,
      std::to_string(omitted_path_count) + " path(s) omitted by ROS output limits",
      response_bytes, configuration.planning.max_response_bytes);
  }
  result_info.returned_path_count = boundedUint32(response->path_results.size());
  result_info.output_complete =
    omitted_path_count == 0 &&
    response->path_results.size() == result.path_solutions.size();
  if (!result_info.output_complete)
    markPathOutputIncomplete(result_info, path_output_limits);
  result_info.request_satisfied =
    result_info.requested_path_count > 0 &&
    result_info.returned_path_count == result_info.requested_path_count;
  response->success = !response->path_results.empty();

  const size_t response_budget_remaining =
    configuration.planning.max_response_bytes > response_bytes ?
    configuration.planning.max_response_bytes - response_bytes : 0;
  const size_t debug_by_bytes = response_budget_remaining / kDebugNodeBytes;
  const size_t debug_node_count = request->include_debug ?
    std::min({
      nodes.size(), configuration.planning.max_debug_nodes, debug_by_bytes}) :
    0U;
  response->debug_nodes.reserve(debug_node_count);
  for (size_t i = 0; i < debug_node_count; ++i)
  {
    if ((i & 0xffu) == 0u && stop_requested && stop_requested())
    {
      initializePlanningResponse(
        *response, request->k, map_id, request->include_debug);
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
  if (request->include_debug &&
      emitted_debug_nodes < result.expanded_nodes)
  {
    appendResponseNotice(response->message,
      "Debug output limited to " + std::to_string(emitted_debug_nodes) +
      " node(s)", response_bytes,
      configuration.planning.max_response_bytes);
  }

  result_info.debug_output_complete =
    !request->include_debug || emitted_debug_nodes == result.expanded_nodes;
  response->success = !response->path_results.empty();

  if (stop_requested && stop_requested())
  {
    initializePlanningResponse(
      *response, request->k, map_id, request->include_debug);
    result_info.found_path_count = boundedUint32(result.path_solutions.size());
    result_info.expanded_nodes = boundedUint64(result.expanded_nodes);
    result_info.map_time_ms = result.map_time_ms;
    result_info.plan_time_ms = result.plan_time_ms;
    markCancelled(*response);
    response->message = "Planning was cancelled before publishing results";
    return;
  }

  if (response->success && result.polymap)
  {
    std::vector<PathSolution> visualization_solutions;
    visualization_solutions.reserve(emitted_solution_indices.size());
    for (const size_t solution_index : emitted_solution_indices)
    {
      visualization_solutions.emplace_back(
        std::move(result.path_solutions[solution_index]));
    }
    last_frame_id_ = grid.header.frame_id;

    publishNonHomotopicPaths(
      visualization_solutions, work_map, grid.header.frame_id,
      configuration.planning.max_path_points,
      configuration.planning.max_response_bytes);
    if (request->include_debug)
    {
      publishPolyObstacles(*result.polymap, work_map, grid.header.frame_id,
        configuration.planning.max_response_bytes);
      publishCDT(*result.polymap, work_map, grid.header.frame_id,
        configuration.planning.max_response_bytes);
      publishDebugTree(
        nodes, work_map, grid.header.frame_id,
        configuration.planning.max_debug_nodes,
        configuration.planning.max_response_bytes);
    }
  }

}
catch (const std::exception& exception)
{
  std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
  clearVisualizationsLocked();
  initializePlanningResponse(
    response_value, request_value.k, map_id, request_value.include_debug);
  if (stop_requested && stop_requested())
  {
    markCancelled(response_value);
    response_value.message = "Planning was cancelled while handling an exception";
  }
  else
  {
    response_value.result_info.status = PlanningResultInfo::STATUS_FAILED;
    setBoundedExceptionMessage(response_value.message, exception.what());
  }
  RCLCPP_ERROR(get_logger(), "%s", response_value.message.c_str());
}
catch (...)
{
  std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
  clearVisualizationsLocked();
  initializePlanningResponse(
    response_value, request_value.k, map_id, request_value.include_debug);
  if (stop_requested && stop_requested())
  {
    markCancelled(response_value);
    response_value.message = "Planning was cancelled while handling an unknown exception";
  }
  else
  {
    response_value.result_info.status = PlanningResultInfo::STATUS_FAILED;
    response_value.message = "Raystar request failed with an unknown exception";
  }
  RCLCPP_ERROR(get_logger(), "%s", response_value.message.c_str());
}

bool RaystarNode::buildPathMsg(
  const PathSolution& solution, const GridMap& grid_map,
  const std::string& frame_id, size_t max_path_points,
  nav_msgs::msg::Path& msg, std::string& error) const
{
  msg = nav_msgs::msg::Path{};
  error.clear();
  msg.header.stamp = now();
  msg.header.frame_id = frame_id;

  std::vector<Point2d> interpolated;
  if (!interpolateProjectedPath(
      solution, max_path_points, interpolated, error))
  {
    return false;
  }
  msg.poses.reserve(interpolated.size());
  for (const auto& point : interpolated)
  {
    geometry_msgs::msg::PoseStamped pose;
    pose.header = msg.header;
    pose.pose.orientation.w = 1.0;
    double wx, wy;
    if (!continuousGridToWorld(grid_map, point, wx, wy))
    {
      error = "path contains a point outside the finite world transform";
      msg.poses.clear();
      return false;
    }
    pose.pose.position.x = wx;
    pose.pose.position.y = wy;
    msg.poses.push_back(pose);
  }

  return true;
}

void RaystarNode::publishPolyObstacles(
  const Polymap& polymap, const GridMap& grid_map,
  const std::string& frame_id, size_t max_marker_bytes)
{
  const auto stamp = now();
  auto array = makeMarkerSnapshot(frame_id, stamp);
  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = frame_id;
  marker.header.stamp = stamp;
  marker.ns = "polygons";
  marker.id = 0;
  marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
  marker.action = visualization_msgs::msg::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.scale.x = 0.05;
  marker.color.r = 1.0;
  marker.color.a = 1.0;

  const size_t topic_budget = markerTopicBudget(max_marker_bytes);
  const size_t point_limit = markerPointLimit(topic_budget);
  const size_t marker_limit = markerEntryLimit(topic_budget);
  size_t emitted_points = 0;
  size_t emitted_markers = 0;
  bool truncated = false;
  for (const auto& ob : polymap.obstacles())
  {
    if (emitted_markers >= marker_limit)
    {
      truncated = true;
      break;
    }
    marker.id++;
    marker.points.clear();
    marker.colors.clear();

    for (auto it = ob.ordered_vertices_.begin();
      it != ob.ordered_vertices_.end(); ++it)
    {
      auto nxt = std::next(it);
      if (nxt == ob.ordered_vertices_.end())
        nxt = ob.ordered_vertices_.begin();

      if (!canAppendCount(emitted_points, 2, point_limit))
      {
        truncated = true;
        break;
      }

      double wx1, wy1, wx2, wy2;
      mapToWorld(grid_map,
        static_cast<unsigned int>(it->first),
        static_cast<unsigned int>(it->second), wx1, wy1);
      mapToWorld(grid_map,
        static_cast<unsigned int>(nxt->first),
        static_cast<unsigned int>(nxt->second), wx2, wy2);

      geometry_msgs::msg::Point p1, p2;
      p1.x = wx1; p1.y = wy1;
      p2.x = wx2; p2.y = wy2;
      marker.points.push_back(p1);
      marker.points.push_back(p2);

      std_msgs::msg::ColorRGBA c;
      c.r = 1.0; c.a = 1.0;
      marker.colors.push_back(c);
      marker.colors.push_back(c);
      emitted_points += 2;
    }
    if (!marker.points.empty())
    {
      array.markers.push_back(marker);
      ++emitted_markers;
    }
    if (truncated)
      break;
  }
  (void)truncated;
  poly_obstacle_pub_->publish(array);
}

void RaystarNode::publishCDT(
  const Polymap& polymap, const GridMap& grid_map,
  const std::string& frame_id, size_t max_marker_bytes)
{
  const size_t topic_budget = markerTopicBudget(max_marker_bytes);
  const size_t point_limit = markerPointLimit(topic_budget);
  const size_t marker_limit = markerEntryLimit(topic_budget);
  const size_t edge_limit = point_limit / 2;
  auto cdt_edges = polymap.getCDTEdges(edge_limit);

  const auto stamp = now();
  auto array = makeMarkerSnapshot(frame_id, stamp);

  // Internal (non-constrained) edges — thin blue lines
  visualization_msgs::msg::Marker int_marker;
  int_marker.header.frame_id = frame_id;
  int_marker.header.stamp = stamp;
  int_marker.ns = "cdt_internal";
  int_marker.id = 0;
  int_marker.type = visualization_msgs::msg::Marker::LINE_LIST;
  int_marker.action = visualization_msgs::msg::Marker::ADD;
  int_marker.pose.orientation.w = 1.0;
  int_marker.scale.x = 0.015;
  int_marker.color.r = 0.3;
  int_marker.color.g = 0.5;
  int_marker.color.b = 0.8;
  int_marker.color.a = 0.4;

  // Constrained (obstacle) edges — thick red lines
  visualization_msgs::msg::Marker con_marker;
  con_marker.header.frame_id = frame_id;
  con_marker.header.stamp = stamp;
  con_marker.ns = "cdt_constrained";
  con_marker.id = 0;
  con_marker.type = visualization_msgs::msg::Marker::LINE_LIST;
  con_marker.action = visualization_msgs::msg::Marker::ADD;
  con_marker.pose.orientation.w = 1.0;
  con_marker.scale.x = 0.06;
  con_marker.color.r = 1.0;
  con_marker.color.g = 0.2;
  con_marker.color.b = 0.2;
  con_marker.color.a = 0.9;

  size_t emitted_points = 0;
  for (const auto& e : cdt_edges)
  {
    if (!canAppendCount(emitted_points, 2, point_limit))
      break;
    geometry_msgs::msg::Point p1, p2;
    // Promote integer vertices before multiplying by the ROS float
    // resolution.  Otherwise the multiplication is performed in float and
    // large grid indices can lose several low bits before being assigned to
    // the double-valued Marker point.
    const double resolution = static_cast<double>(grid_map.resolution);
    p1.x = grid_map.origin_x + static_cast<double>(e.a.first) * resolution;
    p1.y = grid_map.origin_y + static_cast<double>(e.a.second) * resolution;
    p2.x = grid_map.origin_x + static_cast<double>(e.b.first) * resolution;
    p2.y = grid_map.origin_y + static_cast<double>(e.b.second) * resolution;

    if (e.is_constrained)
    {
      con_marker.points.push_back(p1);
      con_marker.points.push_back(p2);
    }
    else
    {
      int_marker.points.push_back(p1);
      int_marker.points.push_back(p2);
    }
    emitted_points += 2;
  }

  if (marker_limit > 0)
    array.markers.push_back(int_marker);
  if (marker_limit > 1)
    array.markers.push_back(con_marker);
  cdt_pub_->publish(array);
}

void RaystarNode::publishNonHomotopicPaths(
  const std::vector<PathSolution>& solutions, const GridMap& grid_map,
  const std::string& frame_id, size_t max_path_points,
  size_t max_marker_bytes)
{
  // This snapshot can be replayed long after planning.  A zero stamp asks
  // RViz for the latest transform instead of pinning the cached message to a
  // time that may have fallen out of a dynamic TF buffer.
  const auto stamp = rclcpp::Time(0, 0, get_clock()->get_clock_type());
  auto array = makeMarkerSnapshot(frame_id, stamp);
  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = frame_id;
  marker.header.stamp = stamp;
  marker.id = 0;
  marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
  marker.action = visualization_msgs::msg::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.scale.x = 0.1;
  marker.color.a = 1.0;

  const size_t topic_budget = markerTopicBudget(max_marker_bytes);
  const size_t point_limit = markerPointLimit(topic_budget);
  const size_t marker_limit = markerEntryLimit(topic_budget);
  size_t emitted_points = 0;
  size_t emitted_markers = 0;
  int num_div = static_cast<int>(std::ceil(std::sqrt(solutions.size())));
  int step = 100 / (num_div + 1);

  for (size_t si = 0; si < solutions.size(); ++si)
  {
    if (emitted_markers >= marker_limit || emitted_points >= point_limit)
      break;
    marker.id++;
    marker.ns = "path_" + std::to_string(si + 1);
    marker.points.clear();
    marker.colors.clear();

    int ri = (si / ((num_div + 1) * (num_div + 1)));
    int gi = (static_cast<int>(si) / (num_div + 1)) % (num_div + 1);
    int bi = static_cast<int>(si) % (num_div + 1);
    std_msgs::msg::ColorRGBA color;
    color.r = (100 + ri * step) / 255.0;
    color.g = (100 + gi * step) / 255.0;
    color.b = (100 + bi * step) / 255.0;
    color.a = 1.0;

    std::vector<Point2d> interpolated;
    std::string interpolation_error;
    const size_t remaining_points = point_limit - emitted_points;
    if (!interpolateProjectedPath(
        solutions[si], std::min(max_path_points, remaining_points),
        interpolated,
        interpolation_error))
    {
      RCLCPP_WARN(get_logger(),
        "Skipping path visualization because output limits rejected it: %s",
        interpolation_error.c_str());
      continue;
    }
    for (const auto& point : interpolated)
    {
      double wx, wy;
      if (!continuousGridToWorld(grid_map, point, wx, wy))
        continue;
      geometry_msgs::msg::Point path_point;
      path_point.x = wx;
      path_point.y = wy;
      marker.points.push_back(path_point);
      marker.colors.push_back(color);
      ++emitted_points;
    }
    if (!marker.points.empty())
    {
      array.markers.push_back(marker);
      ++emitted_markers;
    }
  }
  auto snapshot = std::make_shared<MarkerArray>(std::move(array));
  non_homotopic_pub_->publish(*snapshot);
  if (path_visualization_timer_)
    cached_path_visualization_ = std::move(snapshot);
}

void RaystarNode::clearVisualizationsLocked() noexcept
{
  std::string frame_id;
  try
  {
    // Accepted request frames are bounded by kMaxFrameIdBytes.  If even this
    // small copy cannot be made, an empty frame is still safe for DELETEALL.
    frame_id = last_frame_id_;
  }
  catch (...)
  {
    frame_id.clear();
  }

  // Drop the previous search tree before the next request's map conversion.
  // clear() alone retains each Node/vector capacity and can leave a large
  // failed or invalid request resident indefinitely.
  core_.resetSearchState();
  cached_path_visualization_.reset();
  last_frame_id_.clear();

  // Publishing is normally non-throwing, but a middleware allocation failure
  // must not escape from an exception handler and terminate the executor.  A
  // failure on one topic must also not prevent best-effort clears on the
  // remaining topics.
  try
  {
    const auto clear_snapshot = makeMarkerSnapshot(frame_id, now());
    for (const auto& publisher : {
        non_homotopic_pub_, poly_obstacle_pub_, debug_tree_pub_, cdt_pub_})
    {
      try
      {
        publisher->publish(clear_snapshot);
      }
      catch (...)
      {
        // Continue clearing the other durable topics.  The local caches and
        // Core state have already been released.
      }
    }
  }
  catch (...)
  {
    // Even the small DELETEALL snapshot could not be allocated.
  }
}

void RaystarNode::republishCachedPathVisualization()
{
  std::unique_lock<std::mutex> planner_lock(
    planner_cache_mutex_, std::try_to_lock);
  if (!planner_lock.owns_lock() ||
      planning_busy_.load(std::memory_order_acquire))
  {
    return;
  }
  if (!cached_path_visualization_)
    return;

  // Keep publish inside planner_cache_mutex_.  Copying the shared_ptr and
  // publishing after unlock would allow a new request to publish DELETEALL
  // first and this old snapshot second, resurrecting stale paths.
  try
  {
    if (non_homotopic_pub_->get_subscription_count() == 0)
      return;
    non_homotopic_pub_->publish(*cached_path_visualization_);
  }
  catch (const std::exception& exception)
  {
    RCLCPP_WARN_THROTTLE(get_logger(), *get_clock(), 5000,
      "Could not republish cached path visualization: %s", exception.what());
  }
  catch (...)
  {
    RCLCPP_WARN_THROTTLE(get_logger(), *get_clock(), 5000,
      "Could not republish cached path visualization");
  }
}

void RaystarNode::publishDebugTree(
  const std::vector<raystar::Node>& nodes, const GridMap& grid_map,
  const std::string& frame_id, size_t max_debug_nodes,
  size_t max_marker_bytes)
{
  const auto stamp = now();
  auto array = makeMarkerSnapshot(frame_id, stamp);
  const size_t node_count = std::min(nodes.size(), max_debug_nodes);
  if (node_count == 0)
  {
    debug_tree_pub_->publish(array);
    return;
  }

  const size_t topic_budget = markerTopicBudget(max_marker_bytes);
  const size_t point_limit = markerPointLimit(topic_budget);
  const size_t marker_limit = markerEntryLimit(topic_budget);
  // Keep room for the aggregate edge and seed markers appended after the
  // per-node text/region markers.
  const size_t node_marker_limit = marker_limit > 2 ? marker_limit - 2 : 0;
  size_t emitted_points = 0;
  size_t emitted_markers = 0;

  double min_f = std::numeric_limits<double>::max();
  double max_f = std::numeric_limits<double>::lowest();
  const double resolution = static_cast<double>(grid_map.resolution);
  for (size_t i = 0; i < node_count; ++i)
  {
    const auto& n = nodes[i];
    const double f = (n.gCost() + n.hCost()) * resolution;
    if (f < min_f) min_f = f;
    if (f > max_f) max_f = f;
  }

  visualization_msgs::msg::Marker edge_marker;
  edge_marker.header.frame_id = frame_id;
  edge_marker.header.stamp = stamp;
  edge_marker.ns = "tree_edges";
  edge_marker.id = 0;
  edge_marker.type = visualization_msgs::msg::Marker::LINE_LIST;
  edge_marker.action = visualization_msgs::msg::Marker::ADD;
  edge_marker.pose.orientation.w = 1.0;
  edge_marker.scale.x = 0.02;
  edge_marker.color.r = 0.5;
  edge_marker.color.g = 0.5;
  edge_marker.color.b = 0.5;
  edge_marker.color.a = 0.6;

  visualization_msgs::msg::Marker seed_marker;
  seed_marker.header.frame_id = frame_id;
  seed_marker.header.stamp = stamp;
  seed_marker.ns = "seeds";
  seed_marker.id = 0;
  seed_marker.type = visualization_msgs::msg::Marker::POINTS;
  seed_marker.action = visualization_msgs::msg::Marker::ADD;
  seed_marker.pose.orientation.w = 1.0;
  seed_marker.scale.x = 0.12;
  seed_marker.scale.y = 0.12;

  visualization_msgs::msg::Marker text_marker;
  text_marker.header.frame_id = frame_id;
  text_marker.header.stamp = stamp;
  text_marker.ns = "costs";
  text_marker.type = visualization_msgs::msg::Marker::TEXT_VIEW_FACING;
  text_marker.action = visualization_msgs::msg::Marker::ADD;
  text_marker.pose.orientation.w = 1.0;
  text_marker.scale.z = 0.1;
  text_marker.color.r = 1.0;
  text_marker.color.g = 1.0;
  text_marker.color.b = 1.0;
  text_marker.color.a = 1.0;

  for (size_t i = 0; i < node_count; ++i)
  {
    if (emitted_markers >= node_marker_limit ||
        emitted_points >= point_limit)
      break;
    const auto& n = nodes[i];
    double wx, wy;
    mapToWorld(grid_map,
      static_cast<unsigned int>(n.seed().first),
      static_cast<unsigned int>(n.seed().second), wx, wy);

    const double f = (n.gCost() + n.hCost()) * resolution;
    double t = (max_f > min_f) ? (f - min_f) / (max_f - min_f) : 0.0;

    std_msgs::msg::ColorRGBA color;
    color.r = t;
    color.g = 1.0 - t;
    color.b = 0.2;
    color.a = 1.0;

    geometry_msgs::msg::Point p;
    p.x = wx;
    p.y = wy;
    if (!canAppendCount(emitted_points, 1, point_limit))
      break;
    seed_marker.points.push_back(p);
    seed_marker.colors.push_back(color);
    ++emitted_points;

    if (n.parentIndex() >= 0 &&
        static_cast<size_t>(n.parentIndex()) < node_count)
    {
      const auto& parent = nodes[static_cast<size_t>(n.parentIndex())];
      double pwx, pwy;
      mapToWorld(grid_map,
        static_cast<unsigned int>(parent.seed().first),
        static_cast<unsigned int>(parent.seed().second), pwx, pwy);
      geometry_msgs::msg::Point pp;
      pp.x = pwx;
      pp.y = pwy;
      if (canAppendCount(emitted_points, 2, point_limit))
      {
        edge_marker.points.push_back(pp);
        edge_marker.points.push_back(p);
        emitted_points += 2;
      }
    }

    text_marker.id = static_cast<int>(i);
    text_marker.pose.position.x = wx;
    text_marker.pose.position.y = wy + 0.2;
    text_marker.pose.position.z = 0.0;
    std::ostringstream oss;
    oss << "N" << n.index() << " G=" << std::fixed << std::setprecision(1)
        << n.gCost() * resolution
        << " F=" << std::fixed << std::setprecision(1) << f;
    text_marker.text = oss.str();
    if (emitted_markers >= node_marker_limit)
      break;
    array.markers.push_back(text_marker);
    ++emitted_markers;

    if (!n.visibility().empty() && emitted_markers < node_marker_limit &&
        emitted_points < point_limit)
    {
      visualization_msgs::msg::Marker visreg_marker;
      visreg_marker.header.frame_id = frame_id;
      visreg_marker.header.stamp = stamp;
      visreg_marker.ns = "node_" + std::to_string(n.index());
      visreg_marker.id = 0;
      visreg_marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
      visreg_marker.action = visualization_msgs::msg::Marker::ADD;
      visreg_marker.pose.orientation.w = 1.0;
      visreg_marker.scale.x = 0.03;
      visreg_marker.color.g = 0.8;
      visreg_marker.color.a = 0.5;
      for (const auto& v : n.visibility())
      {
        if (!canAppendCount(emitted_points, 1, point_limit))
          break;
        geometry_msgs::msg::Point vp;
        vp.x = grid_map.origin_x + v.first * grid_map.resolution;
        vp.y = grid_map.origin_y + v.second * grid_map.resolution;
        visreg_marker.points.push_back(vp);
        ++emitted_points;
      }
      if (n.index() > 0)
      {
        if (canAppendCount(emitted_points, 1, point_limit))
        {
          geometry_msgs::msg::Point sp;
          sp.x = grid_map.origin_x + static_cast<double>(n.seed().first) * grid_map.resolution;
          sp.y = grid_map.origin_y + static_cast<double>(n.seed().second) * grid_map.resolution;
          visreg_marker.points.push_back(sp);
          ++emitted_points;
        }
      }
      if (visreg_marker.points.size() > 2 &&
          canAppendCount(emitted_points, 1, point_limit))
      {
        visreg_marker.points.push_back(visreg_marker.points.front());
        ++emitted_points;
      }
      if (!visreg_marker.points.empty() && emitted_markers < node_marker_limit)
      {
        array.markers.push_back(visreg_marker);
        ++emitted_markers;
      }

      if (!n.fullVisibility().empty() &&
          n.fullVisibility().size() != n.visibility().size() &&
          emitted_markers < node_marker_limit && emitted_points < point_limit)
      {
        visualization_msgs::msg::Marker full_visreg_marker;
        full_visreg_marker.header.frame_id = frame_id;
        full_visreg_marker.header.stamp = stamp;
        full_visreg_marker.ns = "full_visreg_" + std::to_string(n.index());
        full_visreg_marker.id = 0;
        full_visreg_marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
        full_visreg_marker.action = visualization_msgs::msg::Marker::ADD;
        full_visreg_marker.pose.orientation.w = 1.0;
        full_visreg_marker.scale.x = 0.02;
        full_visreg_marker.color.r = 0.8;
        full_visreg_marker.color.g = 0.6;
        full_visreg_marker.color.b = 0.2;
        full_visreg_marker.color.a = 0.35;
        for (const auto& v : n.fullVisibility())
        {
          if (!canAppendCount(emitted_points, 1, point_limit))
            break;
          geometry_msgs::msg::Point vp;
          vp.x = grid_map.origin_x + v.first * grid_map.resolution;
          vp.y = grid_map.origin_y + v.second * grid_map.resolution;
          full_visreg_marker.points.push_back(vp);
          ++emitted_points;
        }
        if (n.index() > 0)
        {
          if (canAppendCount(emitted_points, 1, point_limit))
          {
            geometry_msgs::msg::Point sp;
            sp.x = grid_map.origin_x + static_cast<double>(n.seed().first) * grid_map.resolution;
            sp.y = grid_map.origin_y + static_cast<double>(n.seed().second) * grid_map.resolution;
            full_visreg_marker.points.push_back(sp);
            ++emitted_points;
          }
        }
        if (full_visreg_marker.points.size() > 2 &&
            canAppendCount(emitted_points, 1, point_limit))
        {
          full_visreg_marker.points.push_back(full_visreg_marker.points.front());
          ++emitted_points;
        }
        if (!full_visreg_marker.points.empty() &&
            emitted_markers < node_marker_limit)
        {
          array.markers.push_back(full_visreg_marker);
          ++emitted_markers;
        }
      }

      if (n.localShortestPath().size() >= 2 &&
          emitted_markers < node_marker_limit && emitted_points < point_limit)
      {
        visualization_msgs::msg::Marker tpath_marker;
        tpath_marker.header.frame_id = frame_id;
        tpath_marker.header.stamp = stamp;
        tpath_marker.ns = "node_" + std::to_string(n.index());
        tpath_marker.id = 1;
        tpath_marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
        tpath_marker.action = visualization_msgs::msg::Marker::ADD;
        tpath_marker.pose.orientation.w = 1.0;
        tpath_marker.scale.x = 0.04;
        tpath_marker.color.r = 0.0;
        tpath_marker.color.g = 0.5;
        tpath_marker.color.b = 1.0;
        tpath_marker.color.a = 0.8;
        for (const auto& wp : n.localShortestPath())
        {
          if (!canAppendCount(emitted_points, 1, point_limit))
            break;
          double tpx, tpy;
          mapToWorld(grid_map,
            static_cast<unsigned int>(wp.first),
            static_cast<unsigned int>(wp.second), tpx, tpy);
          geometry_msgs::msg::Point tp;
          tp.x = tpx;
          tp.y = tpy;
          tpath_marker.points.push_back(tp);
          ++emitted_points;
        }
        if (!tpath_marker.points.empty() && emitted_markers < node_marker_limit)
        {
          array.markers.push_back(tpath_marker);
          ++emitted_markers;
        }
      }
    }
  }

  if (emitted_markers < marker_limit)
  {
    array.markers.push_back(edge_marker);
    ++emitted_markers;
  }
  if (emitted_markers < marker_limit)
    array.markers.push_back(seed_marker);

  debug_tree_pub_->publish(array);
}

}  // namespace raystar
