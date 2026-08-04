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

namespace raystar {

namespace {

// This tolerance validates quaternion encoding only. Rotation components are
// still required to be exactly zero, so no geometric yaw/tilt is ignored.
constexpr double kQuaternionNormTolerance = 1e-12;
constexpr int64_t kDefaultOccupiedThreshold = 99;
constexpr int64_t kDefaultMaxK = 100;
constexpr int64_t kDefaultMaxCostBoundedPaths =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxCostBoundedPaths);
constexpr int64_t kDefaultMaxMultiGoalCount =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxMultiGoalCount);
constexpr int64_t kDefaultMaxTransitionConfigurations =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxTransitionConfigurations);
constexpr int64_t kDefaultMaxTransitionPairs =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxTransitionPairs);
constexpr int64_t kDefaultMaxNodes = 10000;
constexpr int64_t kDefaultPlanningTimeoutMs = 5000;
constexpr int64_t kDefaultMaxMapCells = static_cast<int64_t>(PlanningLimits::kDefaultMaxMapCells);
constexpr int64_t kDefaultMaxMapBytes = static_cast<int64_t>(PlanningLimits::kDefaultMaxMapBytes);
constexpr int64_t kDefaultMaxPathPoints =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxPathPoints);
constexpr int64_t kDefaultMaxDebugNodes =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxDebugNodes);
constexpr int64_t kDefaultMaxResponseBytes =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxResponseBytes);
constexpr int64_t kDefaultPathVisualizationRepublishPeriodMs = 2000;
constexpr int64_t kMaxTimerPeriodMs =
  std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::nanoseconds::max() -
                                                        std::chrono::milliseconds(1))
    .count();
constexpr int64_t kMaxIntParameterValue = static_cast<int64_t>(std::numeric_limits<int>::max());
constexpr int64_t kMaxPlanningTimeoutMs = std::chrono::milliseconds::max().count() - 1;

constexpr int64_t maxSizeTParameterValue() {
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
constexpr size_t kGoalResultBaseBytes = 256;
constexpr size_t kMaxFrameIdBytes = 256;
constexpr size_t kMaxDiagnosticBytes = 768;
constexpr double kPublishedLengthAbsoluteToleranceMeters = 1.0e-9;
constexpr double kPublishedLengthRelativeTolerance = 1.0e-8;

struct IntegerParameterSpec {
  const char* name;
  int64_t default_value;
  int64_t minimum;
  int64_t maximum;
  const char* description;
  bool read_only;
};

constexpr std::array<IntegerParameterSpec, 14> kIntegerParameterSpecs{{
  {"occupied_threshold",
   kDefaultOccupiedThreshold,
   1,
   100,
   "Startup-only occupancy policy: known OccupancyGrid values at or above this threshold are "
   "occupied.",
   true},
  {"max_k",
   kDefaultMaxK,
   1,
   kMaxIntParameterValue,
   "Largest requested number of paths accepted by the server.",
   false},
  {"max_cost_bounded_paths",
   kDefaultMaxCostBoundedPaths,
   1,
   kMaxIntParameterValue,
   "Maximum paths retained while exhaustively enumerating within a cost bound.",
   false},
  {"max_multi_goal_count",
   kDefaultMaxMultiGoalCount,
   1,
   kMaxIntParameterValue,
   "Maximum goals accepted by one shared-tree bounded request.",
   false},
  {"max_transition_configurations",
   kDefaultMaxTransitionConfigurations,
   1,
   kMaxIntParameterValue,
   "Maximum flattened tether configurations accepted by one UPS batch.",
   false},
  {"max_transition_pairs",
   kDefaultMaxTransitionPairs,
   1,
   kMaxIntParameterValue,
   "Maximum directed configuration pairs accepted by one UPS batch.",
   false},
  {"max_nodes",
   kDefaultMaxNodes,
   1,
   kMaxIntParameterValue,
   "Maximum fully constructed search nodes, including the root.",
   false},
  {"planning_timeout_ms",
   kDefaultPlanningTimeoutMs,
   1,
   kMaxPlanningTimeoutMs,
   "Cooperative planning deadline in milliseconds.",
   false},
  {"max_map_cells",
   kDefaultMaxMapCells,
   1,
   kMaxIntParameterValue,
   "Maximum OccupancyGrid width multiplied by height.",
   false},
  {"max_map_bytes",
   kDefaultMaxMapBytes,
   static_cast<int64_t>(kEstimatedPlannerMapBytesPerCell),
   maxSizeTParameterValue(),
   "Conservative map working-set admission budget in bytes.",
   false},
  {"max_path_points",
   kDefaultMaxPathPoints,
   2,
   kMaxIntParameterValue,
   "Maximum total sampled path poses in one response.",
   false},
  {"max_debug_nodes",
   kDefaultMaxDebugNodes,
   0,
   kMaxIntParameterValue,
   "Maximum structured debug nodes; zero disables debug-node output.",
   false},
  {"max_response_bytes",
   kDefaultMaxResponseBytes,
   static_cast<int64_t>(kMinimumResponseBytes),
   maxSizeTParameterValue(),
   "Conservative response and visualization payload budget in bytes.",
   false},
  {"path_visualization_republish_period_ms",
   kDefaultPathVisualizationRepublishPeriodMs,
   0,
   kMaxTimerPeriodMs,
   "Startup-only cached path MarkerArray republish period; zero disables it.",
   true},
}};

const IntegerParameterSpec* findIntegerParameterSpec(const std::string& name) {
  const auto found =
    std::find_if(kIntegerParameterSpecs.begin(),
                 kIntegerParameterSpecs.end(),
                 [&name](const IntegerParameterSpec& spec) { return name == spec.name; });
  return found == kIntegerParameterSpecs.end() ? nullptr : &*found;
}

rcl_interfaces::msg::ParameterDescriptor makeParameterDescriptor(const IntegerParameterSpec& spec) {
  rcl_interfaces::msg::IntegerRange range;
  range.from_value = spec.minimum;
  range.to_value = spec.maximum;
  range.step = 1;

  rcl_interfaces::msg::ParameterDescriptor descriptor;
  descriptor.description = spec.description;
  descriptor.read_only = spec.read_only;
  descriptor.integer_range.emplace_back(range);
  descriptor.additional_constraints =
    spec.read_only
      ? "May be overridden at startup but cannot be changed while the node runs."
      : "An accepted update affects only requests that take a later configuration snapshot.";
  return descriptor;
}

rcl_interfaces::msg::SetParametersResult validateParameterChanges(
  const std::vector<rclcpp::Parameter>& parameters) {
  rcl_interfaces::msg::SetParametersResult result;
  result.successful = true;
  for (const auto& parameter : parameters) {
    const IntegerParameterSpec* spec = findIntegerParameterSpec(parameter.get_name());
    // Other components may register their own parameters on the node.  Do not
    // reject values outside Raystar's owned parameter set.
    if (spec == nullptr)
      continue;
    if (parameter.get_type() != rclcpp::ParameterType::PARAMETER_INTEGER) {
      result.successful = false;
      result.reason = "parameter '" + parameter.get_name() + "' must be an integer";
      return result;
    }
    const int64_t value = parameter.as_int();
    if (value < spec->minimum || value > spec->maximum) {
      result.successful = false;
      result.reason = "parameter '" + parameter.get_name() + "' must be between " +
                      std::to_string(spec->minimum) + " and " + std::to_string(spec->maximum);
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

size_t markerTopicBudget(size_t total_budget) {
  return std::max<size_t>(1, total_budget / kVisualizationTopicCount);
}

size_t markerPointLimit(size_t topic_budget) {
  return topic_budget / kMarkerPointBytes;
}

size_t markerEntryLimit(size_t topic_budget) {
  return topic_budget / kMarkerEntryBytes;
}

bool canAppendCount(size_t current, size_t addition, size_t limit) {
  return addition <= limit && current <= limit - addition;
}

visualization_msgs::msg::MarkerArray makeMarkerSnapshot(const std::string& frame_id,
                                                        const rclcpp::Time& stamp) {
  visualization_msgs::msg::MarkerArray array;
  visualization_msgs::msg::Marker clear_marker;
  clear_marker.header.frame_id = frame_id;
  clear_marker.header.stamp = stamp;
  clear_marker.action = visualization_msgs::msg::Marker::DELETEALL;
  array.markers.emplace_back(std::move(clear_marker));
  return array;
}

bool validatePlanarPose(const geometry_msgs::msg::PoseStamped& pose,
                        const std::string& label,
                        const std::string& map_frame,
                        std::string& error) {
  if (pose.header.frame_id.empty()) {
    error = label + " frame_id must not be empty";
    return false;
  }
  if (pose.header.frame_id != map_frame) {
    error = label + " frame_id '" + pose.header.frame_id + "' does not match map frame_id '" +
            map_frame + "'";
    return false;
  }

  const auto& position = pose.pose.position;
  if (!std::isfinite(position.x) || !std::isfinite(position.y) || !std::isfinite(position.z)) {
    error = label + " position must contain only finite coordinates";
    return false;
  }
  if (position.z != 0.0) {
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
bool continuousGridToWorld(const GridMap& grid_map, const Point2d& point, double& wx, double& wy) {
  wx = std::numeric_limits<double>::quiet_NaN();
  wy = std::numeric_limits<double>::quiet_NaN();
  if (!std::isfinite(static_cast<double>(grid_map.resolution)) || grid_map.resolution <= 0.0f ||
      !std::isfinite(grid_map.origin_x) || !std::isfinite(grid_map.origin_y) ||
      !std::isfinite(point.first) || !std::isfinite(point.second)) {
    return false;
  }

  const double converted_x =
    std::fma(point.first, static_cast<double>(grid_map.resolution), grid_map.origin_x);
  const double converted_y =
    std::fma(point.second, static_cast<double>(grid_map.resolution), grid_map.origin_y);
  if (!std::isfinite(converted_x) || !std::isfinite(converted_y))
    return false;

  wx = converted_x;
  wy = converted_y;
  return true;
}

// Keep interpolation in continuous grid coordinates.  In particular, do not
// round samples back to cells: doing so would move a sub-cell endpoint and
// could make the published path differ from the path that Core validated.
bool countInterpolatedPathPoints(const PathSolution& solution,
                                 size_t max_points,
                                 size_t& point_count,
                                 std::string& error) {
  point_count = 0;
  error.clear();
  if (max_points < 2) {
    error = "max_path_points must be at least 2";
    return false;
  }
  if (solution.turning_points_.size() > max_points - 2) {
    error = "path has more essential waypoints than max_path_points=" + std::to_string(max_points);
    return false;
  }
  if (solution.turning_points_.size() > std::numeric_limits<size_t>::max() - 2) {
    error = "path has too many essential waypoints to represent";
    return false;
  }

  std::vector<Point2d> projected;
  projected.reserve(solution.turning_points_.size() + 2);
  projected.emplace_back(solution.start_);
  for (const auto& point : solution.turning_points_) {
    projected.emplace_back(static_cast<double>(point.first), static_cast<double>(point.second));
  }
  projected.emplace_back(solution.goal_);

  // Each segment contributes its first point and the final goal is appended
  // once.  Check the complete count before allocating the interpolation
  // vector; no rejected path can trigger an oversized reserve().
  size_t count = 1;
  for (size_t i = 0; i + 1 < projected.size(); ++i) {
    const auto& previous = projected[i];
    const auto& next = projected[i + 1];
    const double distance = std::hypot(previous.first - next.first, previous.second - next.second);
    if (!std::isfinite(distance)) {
      error = "path contains a non-finite segment length";
      return false;
    }
    const double rounded_count = std::max(1.0, std::ceil(distance));
    if (rounded_count > static_cast<double>(max_points) ||
        rounded_count >= static_cast<double>(std::numeric_limits<size_t>::max())) {
      error = "path interpolation exceeds max_path_points=" + std::to_string(max_points);
      return false;
    }
    const size_t segment_count = static_cast<size_t>(rounded_count);
    if (segment_count > max_points - count) {
      error = "path interpolation exceeds max_path_points=" + std::to_string(max_points);
      return false;
    }
    count += segment_count;
  }

  point_count = count;
  return true;
}

bool interpolateProjectedPath(const PathSolution& solution,
                              size_t max_points,
                              std::vector<Point2d>& interpolated,
                              std::string& error) {
  interpolated.clear();
  size_t point_count = 0;
  if (!countInterpolatedPathPoints(solution, max_points, point_count, error)) {
    return false;
  }

  std::vector<Point2d> projected;
  projected.reserve(solution.turning_points_.size() + 2);
  projected.emplace_back(solution.start_);
  for (const auto& point : solution.turning_points_) {
    projected.emplace_back(static_cast<double>(point.first), static_cast<double>(point.second));
  }
  projected.emplace_back(solution.goal_);

  interpolated.reserve(point_count);
  for (size_t i = 0; i + 1 < projected.size(); ++i) {
    const auto& previous = projected[i];
    const auto& next = projected[i + 1];
    const double distance = std::hypot(previous.first - next.first, previous.second - next.second);
    // Map dimensions are bounded by signed int in Core, so a finite segment
    // has a representable sample count on supported platforms.  Clamp the
    // degenerate case to one sample to avoid a zero divisor.
    const size_t count = static_cast<size_t>(std::max(1.0, std::ceil(distance)));
    for (size_t j = 0; j < count; ++j) {
      const double fraction = static_cast<double>(j) / static_cast<double>(count);
      interpolated.emplace_back(previous.first + (next.first - previous.first) * fraction,
                                previous.second + (next.second - previous.second) * fraction);
    }
  }
  // Appending the final point explicitly preserves the continuous goal as the
  // last path pose instead of relying on a rounded interpolation sample.
  interpolated.emplace_back(projected.back());
  return true;
}

bool checkedAdd(size_t& total, size_t addition, size_t limit) {
  if (addition > limit || total > limit - addition)
    return false;
  total += addition;
  return true;
}

bool estimatePathResponseBytes(size_t point_count, size_t frame_id_bytes, size_t& bytes) {
  if (frame_id_bytes > std::numeric_limits<size_t>::max() - kPoseBaseBytes ||
      frame_id_bytes > std::numeric_limits<size_t>::max() - kPathBaseBytes)
    return false;
  const size_t bytes_per_pose = kPoseBaseBytes + frame_id_bytes;
  if (point_count >
      (std::numeric_limits<size_t>::max() - kPathBaseBytes - frame_id_bytes) / bytes_per_pose) {
    return false;
  }
  bytes = kPathBaseBytes + frame_id_bytes + point_count * bytes_per_pose;
  return true;
}

bool publishedPathLength(const nav_msgs::msg::Path& path,
                         double& rounded_length,
                         double& conservative_upper_bound,
                         const StopPredicate& stop_requested = {}) {
  rounded_length = std::numeric_limits<double>::quiet_NaN();
  conservative_upper_bound = std::numeric_limits<double>::quiet_NaN();
  if (path.poses.size() < 2)
    return false;
  long double accumulated = 0.0L;
  ConservativeBinary64PathLength certificate;
  for (size_t index = 1; index < path.poses.size(); ++index) {
    if (stop_requested && stop_requested())
      return false;
    const auto& first = path.poses[index - 1].pose.position;
    const auto& second = path.poses[index].pose.position;
    const double segment = std::hypot(second.x - first.x, second.y - first.y);
    if (!std::isfinite(segment) ||
        !certificate.addSegment(first.x, first.y, second.x, second.y)) {
      return false;
    }
    accumulated += static_cast<long double>(segment);
  }
  if (!std::isfinite(accumulated) ||
      accumulated > static_cast<long double>(std::numeric_limits<double>::max())) {
    return false;
  }
  rounded_length = static_cast<double>(accumulated);
  return certificate.upperBound(conservative_upper_bound, stop_requested);
}

bool publishedLengthsMatch(double first, double second) {
  const double tolerance =
    std::max(kPublishedLengthAbsoluteToleranceMeters,
             kPublishedLengthRelativeTolerance * std::max(std::abs(first), std::abs(second)));
  return std::abs(first - second) <= tolerance;
}

enum class MetricBoundEligibility { within_bound, outside_bound, invalid };

MetricBoundEligibility classifyPathViewMetricBound(const BoundedPathView& path,
                                                    const GridMap& map,
                                                    double inclusive_bound,
                                                    const StopToken& stop_token) {
  if (!std::isfinite(path.path_cost) || path.path_cost < 0.0 ||
      !std::isfinite(inclusive_bound) || inclusive_bound < 0.0) {
    return MetricBoundEligibility::invalid;
  }

  ConservativeBinary64PathLength certificate;
  double previous_x = std::numeric_limits<double>::quiet_NaN();
  double previous_y = std::numeric_limits<double>::quiet_NaN();
  bool have_previous = false;
  const auto add_point = [&](const Point2d& grid_point) {
    if (stop_token.poll())
      return false;
    double world_x = std::numeric_limits<double>::quiet_NaN();
    double world_y = std::numeric_limits<double>::quiet_NaN();
    if (!continuousGridToWorld(map, grid_point, world_x, world_y))
      return false;
    if (have_previous &&
        !certificate.addSegment(previous_x, previous_y, world_x, world_y)) {
      return false;
    }
    previous_x = world_x;
    previous_y = world_y;
    have_previous = true;
    return true;
  };

  if (!add_point(path.start))
    return MetricBoundEligibility::invalid;
  for (const auto& turning_point : path.turning_points) {
    if (!add_point(
          {static_cast<double>(turning_point.first),
           static_cast<double>(turning_point.second)})) {
      return MetricBoundEligibility::invalid;
    }
  }
  if (!add_point(path.goal))
    return MetricBoundEligibility::invalid;

  double upper_bound = std::numeric_limits<double>::quiet_NaN();
  const auto certificate_stop = [&]() { return stop_token.poll(); };
  if (!certificate.upperBound(upper_bound, certificate_stop))
    return MetricBoundEligibility::invalid;
  const double core_metric_cost = gridCostToMetric(path.path_cost, map.resolution);
  if (!std::isfinite(core_metric_cost) ||
      !publishedLengthsMatch(core_metric_cost, upper_bound)) {
    return MetricBoundEligibility::invalid;
  }
  return upper_bound <= inclusive_bound ? MetricBoundEligibility::within_bound
                                        : MetricBoundEligibility::outside_bound;
}

MetricBoundEligibility classifyTopologyMetricBound(const nav_msgs::msg::Path& topology_path,
                                                    double core_cost,
                                                    double inclusive_bound,
                                                    const StopPredicate& stop_requested = {}) {
  if (!std::isfinite(core_cost) || core_cost < 0.0 || !std::isfinite(inclusive_bound) ||
      inclusive_bound < 0.0) {
    return MetricBoundEligibility::invalid;
  }
  double rounded_length = std::numeric_limits<double>::quiet_NaN();
  double upper_bound = std::numeric_limits<double>::quiet_NaN();
  if (!publishedPathLength(
        topology_path, rounded_length, upper_bound, stop_requested) ||
      !publishedLengthsMatch(core_cost, upper_bound)) {
    return MetricBoundEligibility::invalid;
  }
  return upper_bound <= inclusive_bound ? MetricBoundEligibility::within_bound
                                        : MetricBoundEligibility::outside_bound;
}

bool finalizePublishedPathResult(
  raystar_interfaces::msg::PathResult& result,
  double core_cost,
  double inclusive_bound = std::numeric_limits<double>::infinity(),
  const StopPredicate& stop_requested = {}) {
  if (!std::isfinite(core_cost) || core_cost < 0.0 || std::isnan(inclusive_bound))
    return false;
  double dense_length = std::numeric_limits<double>::quiet_NaN();
  double topology_length = std::numeric_limits<double>::quiet_NaN();
  double dense_upper_bound = std::numeric_limits<double>::quiet_NaN();
  double topology_upper_bound = std::numeric_limits<double>::quiet_NaN();
  if (!publishedPathLength(
        result.path, dense_length, dense_upper_bound, stop_requested) ||
      !publishedPathLength(result.topology_path,
                           topology_length,
                           topology_upper_bound,
                           stop_requested)) {
    return false;
  }

  // `cost` certifies the geometry actually serialized on the wire. Core's
  // grid-space ceiling can acquire one additional rounding ULP when scaled by
  // a non-dyadic float resolution, even though both serialized world
  // polylines fit the metric bound. Keep it as a faithfulness cross-check, but
  // do not let that representation artifact reject an exact world-space path.
  result.cost = std::max(dense_upper_bound, topology_upper_bound);
  return std::isfinite(result.cost) && result.cost <= inclusive_bound &&
         dense_upper_bound <= inclusive_bound && topology_upper_bound <= inclusive_bound &&
         publishedLengthsMatch(core_cost, result.cost) &&
         publishedLengthsMatch(dense_length, result.cost) &&
         publishedLengthsMatch(topology_length, result.cost);
}

bool finalizeCostBoundedPublishedPathResult(raystar_interfaces::msg::PathResult& result,
                                            double core_cost,
                                            double inclusive_bound,
                                            const StopPredicate& stop_requested = {}) {
  if (finalizePublishedPathResult(
        result, core_cost, inclusive_bound, stop_requested))
    return true;

  // Interpolation is visualization-only. Rounding its intermediate world
  // coordinates can make the serialized dense polyline microscopically
  // longer than the exact topology polyline. At an inclusive one-ULP bound,
  // retry with the collision-free topology poses themselves instead of
  // discarding a mathematically eligible path or weakening the certificate.
  result.path = result.topology_path;
  return finalizePublishedPathResult(result, core_cost, inclusive_bound, stop_requested);
}

bool hasMetricFaithfulWorldTransform(double origin_x,
                                    double origin_y,
                                    double extent_x,
                                    double extent_y,
                                    double resolution) {
  // The published path contract allows a small representation discrepancy
  // from Core's metric cost.  Admit only transforms whose binary64 quantum is
  // comfortably below that budget at both the multiply and add magnitudes.
  // Per-path validation below remains the authoritative final check.
  const double representation_budget =
    std::max(kPublishedLengthAbsoluteToleranceMeters,
             kPublishedLengthRelativeTolerance * resolution) /
    32.0;
  // The representation tolerance alone is insufficient for very small map
  // resolutions: a large origin can have an ULP smaller than 1e-9 metres yet
  // still collapse many adjacent cells. Requiring the coordinate quantum to
  // be much smaller than one cell makes every rounded adjacent transform
  // strictly monotone as well as metrically distinguishable.
  const double coordinate_budget =
    std::min(representation_budget, resolution / 32.0);
  const double max_world_x = std::fma(1.0, extent_x, origin_x);
  const double max_world_y = std::fma(1.0, extent_y, origin_y);
  const double first_world_x = std::fma(1.0, resolution, origin_x);
  const double first_world_y = std::fma(1.0, resolution, origin_y);
  const std::array<double, 8> transform_values{origin_x,
                                               origin_y,
                                               extent_x,
                                               extent_y,
                                               first_world_x,
                                               first_world_y,
                                               max_world_x,
                                               max_world_y};
  if (!std::isfinite(coordinate_budget) || coordinate_budget <= 0.0 ||
      !std::all_of(transform_values.begin(), transform_values.end(), [&](double value) {
        return std::isfinite(value) &&
               binary64SpacingForSearchPadding(value) <= coordinate_budget;
      })) {
    return false;
  }
  if (!(first_world_x > origin_x) || !(first_world_y > origin_y) ||
      !(max_world_x > origin_x) || !(max_world_y > origin_y)) {
    return false;
  }
  const auto distance_is_represented = [](double actual, double expected) {
    if (!std::isfinite(actual) || actual <= 0.0)
      return false;
    const long double error =
      std::abs(static_cast<long double>(actual) - static_cast<long double>(expected));
    const long double tolerance = std::max(
      static_cast<long double>(kPublishedLengthAbsoluteToleranceMeters),
      static_cast<long double>(kPublishedLengthRelativeTolerance) *
        static_cast<long double>(expected));
    return error <= tolerance;
  };
  return distance_is_represented(first_world_x - origin_x, resolution) &&
         distance_is_represented(first_world_y - origin_y, resolution) &&
         distance_is_represented(max_world_x - origin_x, extent_x) &&
         distance_is_represented(max_world_y - origin_y, extent_y);
}

bool appendResponseNotice(std::string& message,
                          const std::string& notice,
                          size_t& response_bytes,
                          size_t max_response_bytes) {
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

void setBoundedExceptionMessage(std::string& message, const char* detail) {
  constexpr const char* prefix = "Raystar request failed with exception: ";
  message.assign(prefix);
  if (message.size() >= kMaxDiagnosticBytes || detail == nullptr) {
    message.resize(std::min(message.size(), kMaxDiagnosticBytes));
    return;
  }

  size_t detail_bytes = 0;
  const size_t remaining = kMaxDiagnosticBytes - message.size();
  while (detail_bytes < remaining && detail[detail_bytes] != '\0') ++detail_bytes;
  message.append(detail, detail_bytes);
}

uint32_t boundedUint32(size_t value) {
  return static_cast<uint32_t>(
    std::min(value, static_cast<size_t>(std::numeric_limits<uint32_t>::max())));
}

uint64_t boundedUint64(size_t value) {
  if constexpr (sizeof(size_t) > sizeof(uint64_t)) {
    return static_cast<uint64_t>(
      std::min(value, static_cast<size_t>(std::numeric_limits<uint64_t>::max())));
  }
  return static_cast<uint64_t>(value);
}

uint16_t planningLimitMask(PlanningLimitReached limit) {
  switch (limit) {
  case PlanningLimitReached::none:
    return PlanningResultInfo::LIMIT_NONE;
  case PlanningLimitReached::timeout:
    return PlanningResultInfo::LIMIT_TIMEOUT;
  case PlanningLimitReached::max_nodes:
    return PlanningResultInfo::LIMIT_MAX_NODES;
  case PlanningLimitReached::max_path_points:
    return PlanningResultInfo::LIMIT_MAX_PATH_POINTS;
  case PlanningLimitReached::max_paths:
    return PlanningResultInfo::LIMIT_MAX_PATHS;
  case PlanningLimitReached::cancelled:
    return PlanningResultInfo::LIMIT_CANCELLED;
  }
  return PlanningResultInfo::LIMIT_NONE;
}

void markPathOutputIncomplete(PlanningResultInfo& info,
                              uint16_t limit = PlanningResultInfo::LIMIT_NONE) {
  info.output_complete = false;
  info.limits_reached = static_cast<uint16_t>(info.limits_reached | limit);
  // Search interruption and terminal/error outcomes take precedence in the
  // high-level status. The independent output flag and limit bit still retain
  // an additional ROS serialization truncation.
  if (info.status == PlanningResultInfo::STATUS_COMPLETE ||
      info.status == PlanningResultInfo::STATUS_FEWER_PATHS ||
      info.status == PlanningResultInfo::STATUS_NO_PATH) {
    info.status = PlanningResultInfo::STATUS_PARTIAL_OUTPUT;
  }
}

template <typename ResponseT>
void resetPlanningResponse(ResponseT& response) {
  response.success = false;
  response.result_info = PlanningResultInfo{};
  response.message.clear();
  response.path_results.clear();
  response.debug_nodes.clear();
}

template <typename ResponseT>
void setRequestedPathCount(ResponseT& response, int32_t requested) {
  response.result_info.requested_path_count = requested > 0 ? static_cast<uint32_t>(requested) : 0u;
}

template <typename ResponseT>
void initializePlanningResponse(ResponseT& response,
                                uint8_t search_mode,
                                int32_t requested,
                                double requested_max_path_length,
                                const raystar_interfaces::MapId& map_id,
                                bool debug_requested) {
  resetPlanningResponse(response);
  response.result_info.map_id = map_id;
  response.result_info.debug_requested = debug_requested;
  response.result_info.debug_output_complete = !debug_requested;
  response.result_info.search_mode = search_mode;
  response.result_info.requested_max_path_length = requested_max_path_length;
  setRequestedPathCount(response, requested);
}

template <typename ResponseT>
void resetPlanningResponsePreservingRequestMetadata(ResponseT& response) {
  const auto map_id = response.result_info.map_id;
  const auto environment_id = response.result_info.environment_id;
  const uint8_t search_mode = response.result_info.search_mode;
  const uint32_t requested_path_count = response.result_info.requested_path_count;
  const double requested_max_path_length = response.result_info.requested_max_path_length;
  const bool debug_requested = response.result_info.debug_requested;
  resetPlanningResponse(response);
  response.result_info.map_id = map_id;
  response.result_info.environment_id = environment_id;
  response.result_info.search_mode = search_mode;
  response.result_info.requested_path_count = requested_path_count;
  response.result_info.requested_max_path_length = requested_max_path_length;
  response.result_info.debug_requested = debug_requested;
  response.result_info.debug_output_complete = !debug_requested;
}

template <typename ResponseT>
void markCancelled(ResponseT& response) {
  response.success = false;
  response.path_results.clear();
  response.debug_nodes.clear();
  response.result_info.status = PlanningResultInfo::STATUS_CANCELLED;
  response.result_info.limits_reached = static_cast<uint16_t>(response.result_info.limits_reached |
                                                              PlanningResultInfo::LIMIT_CANCELLED);
  response.result_info.request_satisfied = false;
  response.result_info.search_complete = false;
  response.result_info.output_complete = false;
  response.result_info.debug_output_complete = false;
  response.result_info.returned_path_count = 0;
}

template <typename ResponseT>
void markFailed(ResponseT& response) {
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

class PlanningSlotRelease {
public:
  explicit PlanningSlotRelease(std::atomic<bool>& busy) : busy_(busy) {}
  PlanningSlotRelease(const PlanningSlotRelease&) = delete;
  PlanningSlotRelease& operator=(const PlanningSlotRelease&) = delete;
  ~PlanningSlotRelease() {
    busy_.store(false, std::memory_order_release);
  }

private:
  std::atomic<bool>& busy_;
};

class SearchStateRelease {
public:
  explicit SearchStateRelease(RaystarCore& core) : core_(core) {}
  SearchStateRelease(const SearchStateRelease&) = delete;
  SearchStateRelease& operator=(const SearchStateRelease&) = delete;
  ~SearchStateRelease() {
    core_.resetSearchState();
  }

private:
  RaystarCore& core_;
};

template <typename Callback>
class ScopeExit {
public:
  explicit ScopeExit(Callback callback) : callback_(std::move(callback)) {}
  ScopeExit(const ScopeExit&) = delete;
  ScopeExit& operator=(const ScopeExit&) = delete;
  ~ScopeExit() noexcept {
    try {
      callback_();
    } catch (...) {}
  }

private:
  Callback callback_;
};

}  // namespace

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
  configuration.planning.max_cost_bounded_paths =
    static_cast<size_t>(max_cost_bounded_paths);
  configuration.planning.max_multi_goal_count = static_cast<size_t>(max_multi_goal_count);
  configuration.planning.max_transition_configurations =
    static_cast<size_t>(max_transition_configurations);
  configuration.planning.max_transition_pairs =
    static_cast<size_t>(max_transition_pairs);
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
    status.environment_identity_version =
      raystar_interfaces::kEnvironmentIdentityVersion;
    status.occupancy_semantics_version = raystar_interfaces::kOccupancySemanticsVersion;
    status.geometry_semantics_version = raystar_interfaces::kGeometrySemanticsVersion;
    status.topology_semantics_version = raystar_interfaces::kTopologySemanticsVersion;
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
  for (const auto& goal : goals)
    cache_goals.emplace_back(goal);
  return transition_environment_cache_.find(
    TransitionEnvironmentKey(cached_map_.generation,
                             policy,
                             TransitionEnvironmentEndpoint(base),
                             std::move(cache_goals)));
}

void RaystarNode::cacheCompletedTransitionEnvironment(
  const nav_msgs::msg::OccupancyGrid& grid,
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
  for (const auto& goal : goals)
    cache_goals.emplace_back(goal);
  transition_environment_cache_.store(
    TransitionEnvironmentKey(cached_map_.generation,
                             policy,
                             TransitionEnvironmentEndpoint(base),
                             std::move(cache_goals)),
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
  const double extent_x =
    gridCostToMetric(static_cast<double>(width), grid.info.resolution);
  const double extent_y =
    gridCostToMetric(static_cast<double>(height), grid.info.resolution);
  const double max_world_x = std::fma(1.0, extent_x, origin.position.x);
  const double max_world_y = std::fma(1.0, extent_y, origin.position.y);
  if (!std::isfinite(max_world_x) || !std::isfinite(max_world_y)) {
    error = "Invalid map: world-coordinate extent is not finite";
    return false;
  }
  if (!hasMetricFaithfulWorldTransform(origin.position.x,
                                       origin.position.y,
                                       extent_x,
                                       extent_y,
                                       metric_resolution)) {
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

void RaystarNode::handleActionAccepted(const std::shared_ptr<PlanGoalHandle> goal_handle) {
  std::shared_ptr<std::atomic<bool>> cancel_requested;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (goal_handle && active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
      cancel_requested = active_goal_cancel_;
    }
  }

  if (!goal_handle || !cancel_requested) {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    active_goal_reserved_ = false;
    active_goal_cancel_.reset();
    planning_busy_.store(false, std::memory_order_release);
    return;
  }

  bool queued = false;
  {
    std::lock_guard<std::mutex> lock(action_worker_mutex_);
    if (!stop_action_worker_ && !pending_action_job_) {
      pending_action_job_.emplace(ActionJob{goal_handle, cancel_requested});
      queued = true;
    }
  }
  if (queued) {
    action_worker_cv_.notify_one();
    return;
  }

  // Release admission before doing any best-effort result allocation. This
  // branch is only reachable during shutdown or an internal invariant
  // violation; even an allocation failure must not leave the planner busy.
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
      active_goal_reserved_ = false;
      active_goal_cancel_.reset();
    }
    planning_busy_.store(false, std::memory_order_release);
  }
  try {
    auto result = std::make_shared<PlanAction::Result>();
    const auto rejected_goal = goal_handle ? goal_handle->get_goal() : nullptr;
    if (rejected_goal) {
      initializePlanningResponse(
        *result,
        rejected_goal->search_mode,
        rejected_goal->k,
        rejected_goal->max_path_length,
        rejected_goal->map_id,
        rejected_goal->include_debug);
    } else {
      resetPlanningResponse(*result);
    }
    result->result_info.status = PlanningResultInfo::STATUS_FAILED;
    result->message = "Raystar Action worker is unavailable";
    goal_handle->abort(result);
  } catch (...) {}
}

rclcpp_action::GoalResponse RaystarNode::handleGoalSetActionGoal(
  const rclcpp_action::GoalUUID& uuid,
  std::shared_ptr<const GoalSetAction::Goal> goal) {
  if (!goal || shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
    return rclcpp_action::GoalResponse::REJECT;

  std::shared_ptr<std::atomic<bool>> cancel_requested;
  try {
    cancel_requested = std::make_shared<std::atomic<bool>>(false);
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Rejecting goal-set Action: could not allocate cancellation state");
    return rclcpp_action::GoalResponse::REJECT;
  }
  bool expected_idle = false;
  if (!planning_busy_.compare_exchange_strong(expected_idle, true, std::memory_order_acq_rel)) {
    RCLCPP_WARN(get_logger(), "Rejecting goal-set Action because the capacity-one planner is busy");
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

rclcpp_action::CancelResponse RaystarNode::handleGoalSetActionCancel(
  const std::shared_ptr<GoalSetGoalHandle> goal_handle) {
  if (!goal_handle)
    return rclcpp_action::CancelResponse::REJECT;
  std::lock_guard<std::mutex> lock(action_state_mutex_);
  if (!active_goal_reserved_ || active_goal_id_ != goal_handle->get_goal_id() ||
      !active_goal_cancel_) {
    return rclcpp_action::CancelResponse::REJECT;
  }
  active_goal_cancel_->store(true, std::memory_order_release);
  return rclcpp_action::CancelResponse::ACCEPT;
}

void RaystarNode::handleGoalSetActionAccepted(
  const std::shared_ptr<GoalSetGoalHandle> goal_handle) {
  std::shared_ptr<std::atomic<bool>> cancel_requested;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (goal_handle && active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id())
      cancel_requested = active_goal_cancel_;
  }
  if (!goal_handle || !cancel_requested) {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    active_goal_reserved_ = false;
    active_goal_cancel_.reset();
    planning_busy_.store(false, std::memory_order_release);
    return;
  }

  bool queued = false;
  {
    std::lock_guard<std::mutex> lock(action_worker_mutex_);
    if (!stop_action_worker_ && !pending_action_job_) {
      pending_action_job_.emplace(GoalSetActionJob{goal_handle, cancel_requested});
      queued = true;
    }
  }
  if (queued) {
    action_worker_cv_.notify_one();
    return;
  }

  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
      active_goal_reserved_ = false;
      active_goal_cancel_.reset();
    }
    planning_busy_.store(false, std::memory_order_release);
  }
  try {
    auto result = std::make_shared<GoalSetAction::Result>();
    result->result_info.map_id = goal_handle->get_goal()->map_id;
    result->result_info.status = PlanningResultInfo::STATUS_FAILED;
    result->message = "Raystar Action worker is unavailable";
    goal_handle->abort(result);
  } catch (...) {}
}

rclcpp_action::GoalResponse RaystarNode::handleTransitionActionGoal(
  const rclcpp_action::GoalUUID& uuid,
  std::shared_ptr<const TransitionAction::Goal> goal) {
  if (!goal || shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
    return rclcpp_action::GoalResponse::REJECT;

  std::shared_ptr<std::atomic<bool>> cancel_requested;
  try {
    cancel_requested = std::make_shared<std::atomic<bool>>(false);
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Rejecting UPS Action: could not allocate cancellation state");
    return rclcpp_action::GoalResponse::REJECT;
  }
  bool expected_idle = false;
  if (!planning_busy_.compare_exchange_strong(expected_idle, true, std::memory_order_acq_rel)) {
    RCLCPP_WARN(get_logger(), "Rejecting UPS Action because the capacity-one planner is busy");
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

rclcpp_action::CancelResponse RaystarNode::handleTransitionActionCancel(
  const std::shared_ptr<TransitionGoalHandle> goal_handle) {
  if (!goal_handle)
    return rclcpp_action::CancelResponse::REJECT;
  std::lock_guard<std::mutex> lock(action_state_mutex_);
  if (!active_goal_reserved_ || active_goal_id_ != goal_handle->get_goal_id() ||
      !active_goal_cancel_) {
    return rclcpp_action::CancelResponse::REJECT;
  }
  active_goal_cancel_->store(true, std::memory_order_release);
  return rclcpp_action::CancelResponse::ACCEPT;
}

void RaystarNode::handleTransitionActionAccepted(
  const std::shared_ptr<TransitionGoalHandle> goal_handle) {
  std::shared_ptr<std::atomic<bool>> cancel_requested;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (goal_handle && active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id())
      cancel_requested = active_goal_cancel_;
  }
  if (!goal_handle || !cancel_requested) {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    active_goal_reserved_ = false;
    active_goal_cancel_.reset();
    planning_busy_.store(false, std::memory_order_release);
    return;
  }

  bool queued = false;
  {
    std::lock_guard<std::mutex> lock(action_worker_mutex_);
    if (!stop_action_worker_ && !pending_action_job_) {
      pending_action_job_.emplace(TransitionActionJob{goal_handle, cancel_requested});
      queued = true;
    }
  }
  if (queued) {
    action_worker_cv_.notify_one();
    return;
  }

  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
      active_goal_reserved_ = false;
      active_goal_cancel_.reset();
    }
    planning_busy_.store(false, std::memory_order_release);
  }
  try {
    auto result = std::make_shared<TransitionAction::Result>();
    result->map_id = goal_handle->get_goal()->map_id;
    result->status = TransitionAction::Result::STATUS_FAILED;
    result->message = "Raystar Action worker is unavailable";
    goal_handle->abort(result);
  } catch (...) {}
}

void RaystarNode::actionWorkerLoop() noexcept {
  while (true) {
    std::optional<PendingActionJob> job;
    {
      std::unique_lock<std::mutex> lock(action_worker_mutex_);
      action_worker_cv_.wait(
        lock, [this]() { return stop_action_worker_ || pending_action_job_.has_value(); });
      if (stop_action_worker_ && !pending_action_job_)
        return;
      job.emplace(std::move(*pending_action_job_));
      pending_action_job_.reset();
    }
    std::visit(
      [this](auto&& typed_job) {
        using Job = std::decay_t<decltype(typed_job)>;
        if constexpr (std::is_same_v<Job, ActionJob>) {
          executeAction(typed_job.goal_handle, typed_job.cancel_requested);
        } else if constexpr (std::is_same_v<Job, GoalSetActionJob>) {
          executeGoalSetAction(typed_job.goal_handle, typed_job.cancel_requested);
        } else {
          executeTransitionAction(typed_job.goal_handle, typed_job.cancel_requested);
        }
      },
      *job);
  }
}

void RaystarNode::executeAction(
  const std::shared_ptr<PlanGoalHandle> goal_handle,
  const std::shared_ptr<std::atomic<bool>>& cancel_requested) noexcept {
  enum class TerminalState { succeeded, aborted, canceled };

  std::shared_ptr<PlanAction::Result> result;
  TerminalState terminal_state = TerminalState::aborted;
  const auto goal_is_canceling = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_canceling();
    } catch (...) {
      return false;
    }
  };
  const auto goal_is_active = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_active();
    } catch (...) {
      return false;
    }
  };

  try {
    result = std::make_shared<PlanAction::Result>();
    resetPlanningResponse(*result);
    const std::weak_ptr<PlanGoalHandle> weak_goal_handle(goal_handle);
    const StopPredicate stop_requested = [this, weak_goal_handle, cancel_requested]() {
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
    if (goal) {
      nav_msgs::msg::OccupancyGrid::ConstSharedPtr cached_map;
      std::string map_error;
      if (!resolveCachedMap(goal->map_id, cached_map, map_error)) {
        initializePlanningResponse(*result,
                                   goal->search_mode,
                                   goal->k,
                                   goal->max_path_length,
                                   goal->map_id,
                                   goal->include_debug);
        result->result_info.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
        result->message = std::move(map_error);
      } else {
        executePlanning(*goal, *result, *cached_map, goal->map_id, stop_requested);
      }
    } else {
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->message = "Raystar Action goal data is unavailable";
    }

    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok()) {
      markFailed(*result);
      if (result->message.empty())
        result->message = "Planning stopped because the Raystar node is shutting down";
      terminal_state = TerminalState::aborted;
    } else if (goal_is_canceling()) {
      markCancelled(*result);
      if (result->message.empty())
        result->message = "Planning was cancelled";
      terminal_state = TerminalState::canceled;
    } else if (result->success) {
      terminal_state = TerminalState::succeeded;
    } else {
      terminal_state = TerminalState::aborted;
    }
  } catch (const std::exception& exception) {
    try {
      if (!result)
        result = std::make_shared<PlanAction::Result>();
      resetPlanningResponsePreservingRequestMetadata(*result);
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      setBoundedExceptionMessage(result->message, exception.what());
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
    RCLCPP_ERROR(get_logger(), "Raystar Action worker failed: %s", exception.what());
  } catch (...) {
    try {
      if (!result)
        result = std::make_shared<PlanAction::Result>();
      resetPlanningResponsePreservingRequestMetadata(*result);
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->message = "Raystar Action worker failed with an unknown exception";
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
    RCLCPP_ERROR(get_logger(), "Raystar Action worker failed with an unknown exception");
  }

  // Linearize completion against handleActionCancel().  If the cancel
  // callback stored its per-goal flag first, retain the reservation while
  // rclcpp_action performs its internal transition after that callback
  // returns.  Otherwise clear the reservation under the same mutex, causing a
  // later cancel request to be rejected before this normal terminal result.
  bool explicit_cancel = false;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    explicit_cancel = cancel_requested->load(std::memory_order_acquire) &&
                      !shutting_down_.load(std::memory_order_acquire);
    if (!explicit_cancel) {
      if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
        terminal_state = TerminalState::aborted;
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
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
    try {
      std::lock_guard<std::mutex> lock(action_state_mutex_);
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
    } catch (...) {
      // Admission still has to be released if middleware state is already
      // tearing down; the node destructor owns the remaining cleanup.
    }
    planning_busy_.store(false, std::memory_order_release);
  });

  if (explicit_cancel) {
    terminal_state = TerminalState::canceled;
    const auto transition_deadline = std::chrono::steady_clock::now() + std::chrono::seconds(1);
    while (rclcpp::ok() && goal_is_active() && !goal_is_canceling() &&
           std::chrono::steady_clock::now() < transition_deadline) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (!goal_is_canceling()) {
      terminal_state = TerminalState::aborted;
      if (result) {
        resetPlanningResponsePreservingRequestMetadata(*result);
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

  if (terminal_state == TerminalState::canceled) {
    markCancelled(*result);
    if (result->message.empty())
      result->message = "Planning was cancelled";
    std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
    clearVisualizationsLocked();
  }

  try {
    if (terminal_state == TerminalState::canceled)
      goal_handle->canceled(result);
    else if (terminal_state == TerminalState::succeeded)
      goal_handle->succeed(result);
    else
      goal_handle->abort(result);
  } catch (const std::exception& exception) {
    // The bounded cancel-transition fallback can race with the transition at
    // its deadline.  If abort() lost that race, finish with the now-valid
    // canceled transition instead of leaving the goal non-terminal.
    if (terminal_state != TerminalState::canceled && goal_is_canceling()) {
      try {
        markCancelled(*result);
        result->message = "Planning was cancelled";
        {
          std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
          clearVisualizationsLocked();
        }
        goal_handle->canceled(result);
        return;
      } catch (...) {}
    }
    RCLCPP_ERROR(
      get_logger(), "Could not publish Raystar Action terminal result: %s", exception.what());
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Could not publish Raystar Action terminal result");
  }
}

void RaystarNode::executeGoalSetAction(
  const std::shared_ptr<GoalSetGoalHandle> goal_handle,
  const std::shared_ptr<std::atomic<bool>>& cancel_requested) noexcept {
  enum class TerminalState { succeeded, aborted, canceled };
  std::shared_ptr<GoalSetAction::Result> result;
  TerminalState terminal_state = TerminalState::aborted;
  const auto goal_is_canceling = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_canceling();
    } catch (...) {
      return false;
    }
  };
  const auto goal_is_active = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_active();
    } catch (...) {
      return false;
    }
  };
  const auto mark_cancelled = [](GoalSetAction::Result& output) {
    output.success = false;
    output.goal_results.clear();
    output.debug_nodes.clear();
    output.result_info.status = PlanningResultInfo::STATUS_CANCELLED;
    output.result_info.limits_reached = static_cast<uint16_t>(
      output.result_info.limits_reached | PlanningResultInfo::LIMIT_CANCELLED);
    output.result_info.request_satisfied = false;
    output.result_info.search_complete = false;
    output.result_info.output_complete = false;
    output.result_info.debug_output_complete = false;
    output.result_info.returned_goal_count = 0;
    output.result_info.completed_goal_count = 0;
    output.result_info.goals_with_paths = 0;
    output.result_info.found_path_count = 0;
    output.result_info.returned_path_count = 0;
  };

  try {
    result = std::make_shared<GoalSetAction::Result>();
    const std::weak_ptr<GoalSetGoalHandle> weak_goal_handle(goal_handle);
    const StopPredicate stop_requested = [this, weak_goal_handle, cancel_requested]() {
      if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
        return true;
      if (!weak_goal_handle.lock())
        return true;
      return cancel_requested->load(std::memory_order_acquire);
    };
    const auto goal = goal_handle->get_goal();
    if (goal) {
      result->result_info.map_id = goal->map_id;
      result->result_info.requested_goal_count = boundedUint32(goal->goals.size());
      result->result_info.debug_requested = goal->include_debug;
      result->result_info.debug_output_complete = !goal->include_debug;
      nav_msgs::msg::OccupancyGrid::ConstSharedPtr cached_map;
      std::string map_error;
      if (!resolveCachedMap(goal->map_id, cached_map, map_error)) {
        result->result_info.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
        result->message = std::move(map_error);
      } else {
        executeGoalSetPlanning(*goal, *result, *cached_map, goal->map_id, stop_requested);
      }
    } else {
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->message = "Raystar goal-set Action data is unavailable";
    }

    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok()) {
      result->success = false;
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->result_info.request_satisfied = false;
      result->result_info.search_complete = false;
      terminal_state = TerminalState::aborted;
    } else if (goal_is_canceling()) {
      mark_cancelled(*result);
      if (result->message.empty())
        result->message = "Planning was cancelled";
      terminal_state = TerminalState::canceled;
    } else if (result->success) {
      terminal_state = TerminalState::succeeded;
    }
  } catch (const std::exception& exception) {
    try {
      if (!result)
        result = std::make_shared<GoalSetAction::Result>();
      result->success = false;
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->goal_results.clear();
      result->debug_nodes.clear();
      setBoundedExceptionMessage(result->message, exception.what());
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
    RCLCPP_ERROR(get_logger(), "Raystar goal-set Action worker failed: %s", exception.what());
  } catch (...) {
    try {
      if (!result)
        result = std::make_shared<GoalSetAction::Result>();
      result->success = false;
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->goal_results.clear();
      result->debug_nodes.clear();
      result->message = "Raystar goal-set Action worker failed with an unknown exception";
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
  }

  bool explicit_cancel = false;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    explicit_cancel = cancel_requested->load(std::memory_order_acquire) &&
                      !shutting_down_.load(std::memory_order_acquire);
    if (!explicit_cancel) {
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
      planning_busy_.store(false, std::memory_order_release);
    }
  }
  ScopeExit explicit_cancel_release([&]() noexcept {
    if (!explicit_cancel)
      return;
    try {
      std::lock_guard<std::mutex> lock(action_state_mutex_);
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
    } catch (...) {}
    planning_busy_.store(false, std::memory_order_release);
  });

  if (explicit_cancel) {
    terminal_state = TerminalState::canceled;
    const auto transition_deadline = std::chrono::steady_clock::now() + std::chrono::seconds(1);
    while (rclcpp::ok() && goal_is_active() && !goal_is_canceling() &&
           std::chrono::steady_clock::now() < transition_deadline) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (!goal_is_canceling()) {
      terminal_state = TerminalState::aborted;
      if (result) {
        result->success = false;
        result->result_info.status = PlanningResultInfo::STATUS_FAILED;
        result->message = "Cancellation was accepted but the Action state transition did not complete";
      }
    }
  }
  if (!result)
    return;
  if (terminal_state == TerminalState::canceled) {
    mark_cancelled(*result);
    if (result->message.empty())
      result->message = "Planning was cancelled";
    std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
    clearVisualizationsLocked();
  }
  try {
    if (terminal_state == TerminalState::canceled)
      goal_handle->canceled(result);
    else if (terminal_state == TerminalState::succeeded)
      goal_handle->succeed(result);
    else
      goal_handle->abort(result);
  } catch (const std::exception& exception) {
    RCLCPP_ERROR(
      get_logger(), "Could not publish Raystar goal-set Action result: %s", exception.what());
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Could not publish Raystar goal-set Action result");
  }
}

void RaystarNode::executeTransitionAction(
  const std::shared_ptr<TransitionGoalHandle> goal_handle,
  const std::shared_ptr<std::atomic<bool>>& cancel_requested) noexcept {
  enum class TerminalState { succeeded, aborted, canceled };
  std::shared_ptr<TransitionAction::Result> result;
  TerminalState terminal_state = TerminalState::aborted;
  const auto goal_is_canceling = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_canceling();
    } catch (...) {
      return false;
    }
  };
  const auto goal_is_active = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_active();
    } catch (...) {
      return false;
    }
  };
  const auto mark_cancelled = [](TransitionAction::Result& output) {
    output.success = false;
    output.status = TransitionAction::Result::STATUS_CANCELLED;
    if (output.message.empty())
      output.message = "UPS transition construction was cancelled";
  };

  try {
    result = std::make_shared<TransitionAction::Result>();
    const std::weak_ptr<TransitionGoalHandle> weak_goal_handle(goal_handle);
    const StopPredicate stop_requested = [this, weak_goal_handle, cancel_requested]() {
      if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
        return true;
      if (!weak_goal_handle.lock())
        return true;
      return cancel_requested->load(std::memory_order_acquire);
    };
    const TransitionProgressCallback progress_callback =
      [weak_goal_handle](std::uint32_t requested,
                         std::uint32_t completed,
                         const char* stage) noexcept {
        try {
          const auto handle = weak_goal_handle.lock();
          if (!handle || !handle->is_active() || !handle->is_executing())
            return;
          auto feedback = std::make_shared<TransitionAction::Feedback>();
          feedback->requested_transition_count = requested;
          feedback->completed_transition_count = completed;
          feedback->stage = stage == nullptr ? "" : stage;
          handle->publish_feedback(feedback);
        } catch (...) {
          // Feedback is best-effort. A canceled/expired handle, allocation
          // failure, or middleware exception must not abort UPS planning.
        }
      };
    const auto goal = goal_handle->get_goal();
    if (goal) {
      result->map_id = goal->map_id;
      result->requested_transition_count = boundedUint32(goal->transition_pairs.size());
      nav_msgs::msg::OccupancyGrid::ConstSharedPtr cached_map;
      std::string map_error;
      if (!resolveCachedMap(goal->map_id, cached_map, map_error)) {
        progress_callback(result->requested_transition_count,
                          0,
                          "validating transition request");
        result->status = TransitionAction::Result::STATUS_INVALID_REQUEST;
        result->message = std::move(map_error);
      } else {
        executeTransitionPlanning(*goal,
                                  *result,
                                  *cached_map,
                                  goal->map_id,
                                  stop_requested,
                                  progress_callback);
      }
    } else {
      result->status = TransitionAction::Result::STATUS_FAILED;
      result->message = "Raystar UPS Action data is unavailable";
    }

    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok()) {
      result->success = false;
      result->status = TransitionAction::Result::STATUS_FAILED;
      terminal_state = TerminalState::aborted;
    } else if (goal_is_canceling()) {
      mark_cancelled(*result);
      terminal_state = TerminalState::canceled;
    } else if (result->success) {
      terminal_state = TerminalState::succeeded;
    }
  } catch (const std::exception& exception) {
    try {
      if (!result)
        result = std::make_shared<TransitionAction::Result>();
      result->success = false;
      result->status = TransitionAction::Result::STATUS_FAILED;
      result->transitions.clear();
      result->completed_transition_count = 0;
      setBoundedExceptionMessage(result->message, exception.what());
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
    RCLCPP_ERROR(get_logger(), "Raystar UPS Action worker failed: %s", exception.what());
  } catch (...) {
    try {
      if (!result)
        result = std::make_shared<TransitionAction::Result>();
      result->success = false;
      result->status = TransitionAction::Result::STATUS_FAILED;
      result->transitions.clear();
      result->completed_transition_count = 0;
      result->message = "Raystar UPS Action worker failed with an unknown exception";
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
  }

  bool explicit_cancel = false;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    explicit_cancel = cancel_requested->load(std::memory_order_acquire) &&
                      !shutting_down_.load(std::memory_order_acquire);
    if (!explicit_cancel) {
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
      planning_busy_.store(false, std::memory_order_release);
    }
  }
  ScopeExit explicit_cancel_release([&]() noexcept {
    if (!explicit_cancel)
      return;
    try {
      std::lock_guard<std::mutex> lock(action_state_mutex_);
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
    } catch (...) {}
    planning_busy_.store(false, std::memory_order_release);
  });

  if (explicit_cancel) {
    terminal_state = TerminalState::canceled;
    const auto transition_deadline = std::chrono::steady_clock::now() + std::chrono::seconds(1);
    while (rclcpp::ok() && goal_is_active() && !goal_is_canceling() &&
           std::chrono::steady_clock::now() < transition_deadline) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (!goal_is_canceling()) {
      terminal_state = TerminalState::aborted;
      if (result) {
        result->success = false;
        result->status = TransitionAction::Result::STATUS_FAILED;
        result->message =
          "Cancellation was accepted but the Action state transition did not complete";
      }
    }
  }
  if (!result)
    return;
  if (terminal_state == TerminalState::canceled)
    mark_cancelled(*result);
  try {
    if (terminal_state == TerminalState::canceled)
      goal_handle->canceled(result);
    else if (terminal_state == TerminalState::succeeded)
      goal_handle->succeed(result);
    else
      goal_handle->abort(result);
  } catch (const std::exception& exception) {
    RCLCPP_ERROR(get_logger(), "Could not publish Raystar UPS Action result: %s", exception.what());
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Could not publish Raystar UPS Action result");
  }
}

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
    const auto metric_admission =
      [&work_map, metric_bound = request->max_path_length](
        const BoundedPathView& path, const StopToken& stop_token) {
        switch (classifyPathViewMetricBound(
          path, work_map, metric_bound, stop_token)) {
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

  // Core has finished every mutable visibility query before returning and
  // exposes only const shared ownership here. Retain that completed geometry,
  // independently of the search tree, for a following exact-key UPS request.
  if (result.polymap && !result.path_solutions.empty()) {
    cacheCompletedTransitionEnvironment(
      grid,
      map_id,
      request->allow_unknown,
      configuration,
      PolymapEndpoint{start_endpoint.cell_.first,
                      start_endpoint.cell_.second,
                      start_endpoint.position_},
      {PolymapEndpoint{
        goal_endpoint.cell_.first, goal_endpoint.cell_.second, goal_endpoint.position_}},
      result.polymap);
  }

  auto& result_info = response->result_info;
  result_info.found_path_count = boundedUint32(result.path_solutions.size());
  result_info.expanded_nodes = boundedUint64(result.expanded_nodes);
  result_info.map_time_ms = result.map_time_ms;
  result_info.plan_time_ms = result.plan_time_ms;
  result_info.limits_reached = planningLimitMask(result.limit_reached);
  result_info.cost_bound_exhausted =
    cost_bounded_mode &&
    (result.completion == PlanningCompletion::cost_bound_exhausted ||
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
  for (size_t solution_index = 0; solution_index < result.path_solutions.size();
       ++solution_index) {
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
    const double core_metric_cost =
      gridCostToMetric(solution.path_cost_, work_map.resolution);
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
  result_info.found_path_count =
    boundedUint32(metric_eligible_solution_indices.size());

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
    if (!buildTopologyPathMsg(sol,
                              work_map,
                              grid.header.frame_id,
                              path_result.topology_path,
                              path_error)) {
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
      return topology_point_count <=
               configuration.planning.max_path_points - emitted_path_points &&
             serialized_path_points <= configuration.planning.max_path_points -
                                           emitted_path_points - topology_point_count;
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
    const auto bytes_fit = [&](size_t serialized_path_points,
                               size_t& candidate_bytes) {
      candidate_bytes = response_bytes;
      return estimatePathResponseBytes(
               serialized_path_points, grid.header.frame_id.size(), path_bytes) &&
             estimatePathResponseBytes(
               topology_point_count, grid.header.frame_id.size(), topology_path_bytes) &&
             checkedAdd(candidate_bytes,
                        path_bytes,
                        configuration.planning.max_response_bytes) &&
             checkedAdd(candidate_bytes,
                        topology_path_bytes,
                        configuration.planning.max_response_bytes);
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
    } else if (!buildPathMsg(
                 sol,
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
    const double core_metric_cost =
      gridCostToMetric(sol.path_cost_, work_map.resolution);
    const bool finalized =
      cost_bounded_mode
        ? finalizeCostBoundedPublishedPathResult(
            path_result, core_metric_cost, request->max_path_length, stop_requested)
        : finalizePublishedPathResult(path_result,
                                      core_metric_cost,
                                      std::numeric_limits<double>::infinity(),
                                      stop_requested);
    if (!finalized) {
      ++omitted_path_count;
      appendResponseNotice(
        response->message,
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
      path_output_limits = static_cast<uint16_t>(
        path_output_limits | PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
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

  if (!stableSortPublishedPathsWithSources(
        response->path_results, emitted_solution_indices)) {
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
  result_info.output_complete = omitted_path_count == 0 &&
                                response->path_results.size() == metric_found_path_count;
  if (!result_info.output_complete)
    markPathOutputIncomplete(result_info, path_output_limits);
  result_info.request_satisfied =
    cost_bounded_mode
      ? result_info.cost_bound_exhausted && result_info.output_complete
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

void RaystarNode::executeGoalSetPlanning(
  const GoalSetAction::Goal& request,
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
    response.message = "Invalid goal-set request: goals and max_path_lengths must have equal lengths";
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
      response.message = "Invalid goal[" + std::to_string(i) +
                         "]: frame_id must be at most 256 bytes";
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
    const auto metric_admission =
      [&work_map, metric_bound = request.max_path_lengths[i]](
        const BoundedPathView& path, const StopToken& stop_token) {
        switch (classifyPathViewMetricBound(
          path, work_map, metric_bound, stop_token)) {
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

  // Reuse is admitted only when every requested protected endpoint produced a
  // path. This proves every endpoint belonged to the reachable set passed to
  // Polymap; otherwise the Core may have built geometry for only a subset and
  // an exact full-request cache key would overstate its construction inputs.
  const bool complete_environment_endpoint_set =
    result.polymap && result.goal_results.size() == core_goals.size() &&
    std::all_of(result.goal_results.begin(), result.goal_results.end(), [](const auto& goal) {
      return !goal.path_solutions.empty();
    });
  if (complete_environment_endpoint_set) {
    std::vector<PolymapEndpoint> environment_goals;
    environment_goals.reserve(core_goals.size());
    for (const auto& goal : core_goals) {
      environment_goals.push_back({goal.endpoint.cell_.first,
                                   goal.endpoint.cell_.second,
                                   goal.endpoint.position_});
    }
    cacheCompletedTransitionEnvironment(
      grid,
      map_id,
      request.allow_unknown,
      configuration,
      PolymapEndpoint{start_endpoint.cell_.first,
                      start_endpoint.cell_.second,
                      start_endpoint.position_},
      environment_goals,
      result.polymap);
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
    aggregate.limits_reached = static_cast<uint16_t>(
      aggregate.limits_reached | PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
    aggregate.output_complete = false;
    return;
  }
  if (!checkedAdd(
        response_bytes, response.message.size(), configuration.planning.max_response_bytes)) {
    response.message = "Response metadata exceeds max_response_bytes=" +
                       std::to_string(configuration.planning.max_response_bytes);
    aggregate.status = PlanningResultInfo::STATUS_PARTIAL_OUTPUT;
    aggregate.limits_reached = static_cast<uint16_t>(
      aggregate.limits_reached | PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
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
    const bool search_complete =
      core_goal_result.outcome == PlanningOutcome::complete ||
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
        !checkedAdd(response_bytes,
                    goal_metadata_bytes,
                    configuration.planning.max_response_bytes)) {
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
    for (size_t solution_index = 0;
         solution_index < core_goal_result.path_solutions.size();
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
      const double core_metric_cost =
        gridCostToMetric(solution.path_cost_, work_map.resolution);
      const auto eligibility = classifyTopologyMetricBound(
        topology_path,
        core_metric_cost,
        request.max_path_lengths[i],
        stop_requested);
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
      if (!buildTopologyPathMsg(solution,
                                work_map,
                                grid.header.frame_id,
                                path_result.topology_path,
                                path_error)) {
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
        info.limits_reached = static_cast<uint16_t>(
          info.limits_reached | PlanningResultInfo::LIMIT_MAX_PATH_POINTS);
        break;
      }
      size_t path_bytes = 0;
      size_t topology_path_bytes = 0;
      size_t candidate_response_bytes = response_bytes;
      const auto bytes_fit = [&](size_t serialized_path_points,
                                 size_t& candidate_bytes) {
        candidate_bytes = response_bytes;
        return estimatePathResponseBytes(
                 serialized_path_points, grid.header.frame_id.size(), path_bytes) &&
               estimatePathResponseBytes(
                 topology_point_count, grid.header.frame_id.size(), topology_path_bytes) &&
               checkedAdd(candidate_bytes,
                          path_bytes,
                          configuration.planning.max_response_bytes) &&
               checkedAdd(candidate_bytes,
                          topology_path_bytes,
                          configuration.planning.max_response_bytes);
      };
      bool response_bytes_fit = bytes_fit(point_count, candidate_response_bytes);
      if (!response_bytes_fit && !topology_only) {
        topology_only = true;
        point_count = topology_point_count;
        response_bytes_fit = bytes_fit(point_count, candidate_response_bytes);
      }
      if (!response_bytes_fit) {
        omitted = metric_found_path_count - goal_output.path_results.size();
        info.limits_reached = static_cast<uint16_t>(
          info.limits_reached | PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
        break;
      }
      if (topology_only) {
        path_result.path = path_result.topology_path;
      } else if (!buildPathMsg(
                   solution,
                   work_map,
                   grid.header.frame_id,
                   configuration.planning.max_path_points - emitted_path_points -
                     topology_point_count,
                   path_result.path,
                   path_error)) {
        ++omitted;
        continue;
      }
      const double core_metric_cost =
        gridCostToMetric(solution.path_cost_, work_map.resolution);
      if (!finalizeCostBoundedPublishedPathResult(
            path_result,
            core_metric_cost,
            request.max_path_lengths[i],
            stop_requested)) {
        ++omitted;
        continue;
      }
      // As in the single-goal path, admission accounting must follow any
      // dense-to-topology fallback performed by the final certificate.
      if (!bytes_fit(path_result.path.poses.size(), candidate_response_bytes)) {
        ++omitted;
        info.limits_reached = static_cast<uint16_t>(
          info.limits_reached | PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES);
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
    if (!stableSortPublishedPathsWithSources(
          goal_output.path_results, goal_emitted_solution_indices)) {
      throw std::logic_error("published goal path/source index cardinality mismatch");
    }
    for (const size_t solution_index : goal_emitted_solution_indices)
      emitted_solutions.emplace_back(i, solution_index);
    info.returned_path_count = boundedUint32(goal_output.path_results.size());
    info.output_complete = omitted == 0 &&
                           goal_output.path_results.size() == metric_found_path_count;
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

  const bool complete_goal_result_set =
    result.goal_results.size() == request.goals.size() &&
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
  const size_t debug_node_count =
    request.include_debug
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
      !raystar_interfaces::environmentIdsEqual(
        request.expected_environment_id, response.environment_id)) {
    response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
    response.message =
      "expected_environment_id does not match the cached map, occupancy policy, "
      "or Raystar planning semantics";
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
    response.message = "UPS transition construction " + stopped_reason() +
                       " before validation";
    return;
  }
  if (request.tether_configurations.empty()) {
    response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
    response.message = "At least one tether configuration is required";
    return;
  }
  if (request.tether_configurations.size() >
      configuration.planning.max_transition_configurations) {
    response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
    response.message = "Tether configuration count exceeds max_transition_configurations=" +
                       std::to_string(configuration.planning.max_transition_configurations);
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
    if (pair.from_configuration >= request.tether_configurations.size() ||
        pair.to_configuration >= request.tether_configurations.size()) {
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
    response.status = effective_stop_requested()
                        ? stopped_status()
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
    gx = canonicalizeWorldGridCoordinate(wx,
                                         work_map.origin_x,
                                         static_cast<double>(work_map.resolution),
                                         work_map.width,
                                         gx);
    gy = canonicalizeWorldGridCoordinate(wy,
                                         work_map.origin_y,
                                         static_cast<double>(work_map.resolution),
                                         work_map.height,
                                         gy);
    if (!std::isfinite(gx) || !std::isfinite(gy) || gx < 0.0 || gy < 0.0 ||
        gx > static_cast<double>(work_map.width) ||
        gy > static_cast<double>(work_map.height)) {
      return false;
    }
    point = {gx, gy};
    return true;
  };

  std::vector<std::vector<Point2d>> configurations;
  configurations.reserve(request.tether_configurations.size());
  std::vector<PolymapEndpoint> endpoints;
  endpoints.reserve(request.tether_configurations.size());
  size_t input_point_count = 0;
  std::optional<ContinuousGridPoint> base;
  for (size_t configuration_index = 0;
       configuration_index < request.tether_configurations.size();
       ++configuration_index) {
    const auto& input_path = request.tether_configurations[configuration_index];
    if (input_path.header.frame_id != grid.header.frame_id ||
        input_path.header.frame_id.size() > kMaxFrameIdBytes) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Tether configuration " + std::to_string(configuration_index) +
                         " Path frame_id does not match the cached map";
      return;
    }
    if (input_path.poses.empty()) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Tether configuration " + std::to_string(configuration_index) +
                         " is empty";
      return;
    }
    if (!canAppendCount(input_point_count,
                        input_path.poses.size(),
                        configuration.planning.max_path_points)) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Input tether paths exceed max_path_points=" +
                         std::to_string(configuration.planning.max_path_points);
      return;
    }
    input_point_count += input_path.poses.size();
    std::vector<Point2d> path;
    path.reserve(input_path.poses.size());
    for (size_t point_index = 0; point_index < input_path.poses.size(); ++point_index) {
      const auto& pose = input_path.poses[point_index];
      std::string pose_error;
      if (!validatePlanarPose(pose,
                              "Tether configuration waypoint",
                              grid.header.frame_id,
                              pose_error)) {
        response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
        response.message = "Configuration " + std::to_string(configuration_index) +
                           " waypoint " + std::to_string(point_index) + ": " + pose_error;
        return;
      }
      Point2d point;
      if (!waypoint_to_grid(pose.pose.position.x, pose.pose.position.y, point)) {
        response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
        response.message = "Configuration " + std::to_string(configuration_index) +
                           " waypoint " + std::to_string(point_index) +
                           " lies outside the map geometry";
        return;
      }
      if (path.empty() || path.back() != point)
        path.push_back(point);
    }

    ContinuousGridPoint configuration_base;
    const auto& first_pose = input_path.poses.front().pose.position;
    if (!worldToContinuousMap(
          work_map, first_pose.x, first_pose.y, configuration_base)) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Tether configuration " + std::to_string(configuration_index) +
                         " base is not a strict map-interior point";
      return;
    }
    if (!base) {
      base = configuration_base;
    } else if (base->x != configuration_base.x || base->y != configuration_base.y) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "All tether configurations must have one identical base point";
      return;
    }
    ContinuousGridPoint endpoint;
    const auto& last_pose = input_path.poses.back().pose.position;
    if (!worldToContinuousMap(work_map, last_pose.x, last_pose.y, endpoint)) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Tether configuration " + std::to_string(configuration_index) +
                         " endpoint is not a strict map-interior point";
      return;
    }
    const PolymapEndpoint polymap_endpoint{static_cast<int>(endpoint.cell_x),
                                           static_cast<int>(endpoint.cell_y),
                                           {endpoint.x, endpoint.y}};
    const auto duplicate_endpoint = std::find_if(
      endpoints.begin(), endpoints.end(), [&polymap_endpoint](const auto& existing) {
        return existing.cell_x == polymap_endpoint.cell_x &&
               existing.cell_y == polymap_endpoint.cell_y &&
               existing.position == polymap_endpoint.position;
      });
    if (duplicate_endpoint == endpoints.end())
      endpoints.push_back(polymap_endpoint);
    configurations.emplace_back(std::move(path));
  }

  const StopToken stop_token(effective_stop_requested);
  const PolymapEndpoint base_endpoint{static_cast<int>(base->cell_x),
                                      static_cast<int>(base->cell_y),
                                      Point2d{base->x, base->y}};
  std::shared_ptr<const Polymap> polymap_owner = findCachedTransitionEnvironment(
    grid, map_id, request.allow_unknown, configuration, base_endpoint, endpoints);
  if (polymap_owner) {
    RCLCPP_DEBUG(get_logger(),
                 "Reusing the completed free-triangle environment for %zu UPS transition(s)",
                 request.transition_pairs.size());
  } else {
    auto polymap_result = endpoints.size() == 1
                            ? Polymap::create(work_map,
                                             base_endpoint.cell_x,
                                             base_endpoint.cell_y,
                                             endpoints.front().cell_x,
                                             endpoints.front().cell_y,
                                             base_endpoint.position,
                                             endpoints.front().position,
                                             stop_token,
                                             configuration.planning)
                            : Polymap::create(work_map,
                                             base_endpoint.cell_x,
                                             base_endpoint.cell_y,
                                             base_endpoint.position,
                                             endpoints,
                                             stop_token,
                                             configuration.planning);
    if (polymap_result.status == PolymapCreateStatus::stopped) {
      response.status = stopped_status();
      response.message = "UPS transition construction " + stopped_reason() +
                         " while building the map";
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
    cacheCompletedTransitionEnvironment(grid,
                                        map_id,
                                        request.allow_unknown,
                                        configuration,
                                        base_endpoint,
                                        endpoints,
                                        polymap_owner);
  }
  const Polymap& polymap = *polymap_owner;

  // Validate every complete configuration before any pairwise shortening.
  // shortenWithinHomotopy() deliberately removes an exact common prefix from
  // alpha_a and alpha_b. Without this admission pass, two caller-supplied
  // references could share an obstacle-crossing prefix that never appears in
  // the composed alpha_a^{-1} * alpha_b path and would therefore evade the
  // pairwise corridor trace.
  progress.publishStage("validating tether configurations");
  for (size_t configuration_index = 0;
       configuration_index < configurations.size();
       ++configuration_index) {
    if (stop_token.poll()) {
      response.status = stopped_status();
      response.message = "UPS transition construction " + stopped_reason() +
                         " while validating tether configuration " +
                         std::to_string(configuration_index);
      return;
    }
    const auto validation = polymap.shortenPathWithinHomotopy(
      configurations[configuration_index], stop_token);
    if (validation.status == HomotopyShorteningStatus::stopped) {
      response.status = stopped_status();
      response.message = "UPS transition construction " + stopped_reason() +
                         " while validating tether configuration " +
                         std::to_string(configuration_index);
      return;
    }
    if (!validation) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Tether configuration " +
                         std::to_string(configuration_index) +
                         " is not a collision-free reference in the cached map";
      if (!validation.message.empty())
        response.message += ": " + validation.message.substr(0, kMaxDiagnosticBytes);
      response.message.resize(std::min(response.message.size(), kMaxDiagnosticBytes));
      return;
    }
    double input_cost = 0.0;
    for (size_t point_index = 1;
         point_index < configurations[configuration_index].size();
         ++point_index) {
      const auto& previous = configurations[configuration_index][point_index - 1];
      const auto& current = configurations[configuration_index][point_index];
      input_cost += std::hypot(current.first - previous.first,
                               current.second - previous.second);
    }
    const double taut_tolerance =
      1.0e-10 * std::max({1.0, input_cost, validation.path_cost});
    if (!std::isfinite(input_cost) ||
        std::abs(input_cost - validation.path_cost) > taut_tolerance) {
      response.status = TransitionAction::Result::STATUS_INVALID_REQUEST;
      response.message = "Tether configuration " +
                         std::to_string(configuration_index) +
                         " is not a locally shortest (taut) reference";
      return;
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
    const auto shortening = RaystarCore::shortenWithinHomotopy(
      polymap,
      configurations[requested_pair.from_configuration],
      configurations[requested_pair.to_configuration],
      stop_token);
    if (shortening.status == HomotopyShorteningStatus::stopped) {
      response.status = stopped_status();
      response.message = "UPS transition construction " + stopped_reason() +
                         " while shortening pair " +
                         std::to_string(index);
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

    if (!canAppendCount(output_point_count,
                        shortening.path.size(),
                        configuration.planning.max_path_points)) {
      response.status = TransitionAction::Result::STATUS_FAILED;
      response.message = "UPS output exceeds max_path_points=" +
                         std::to_string(configuration.planning.max_path_points);
      return;
    }
    size_t path_bytes = 0;
    const bool triangle_bytes_valid =
      output.triangle_occurrences.size() <=
      std::numeric_limits<size_t>::max() / sizeof(uint32_t);
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
      if (!continuousGridToWorld(
            work_map, point, pose.pose.position.x, pose.pose.position.y)) {
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
          !publishedPathLength(output.path,
                               rounded_world_length,
                               certified_world_length,
                               certificate_stop)) {
        const bool stopped = stop_token.poll();
        response.status = stopped ? stopped_status()
                                  : TransitionAction::Result::STATUS_FAILED;
        response.message = stopped
                             ? "UPS transition length certification " + stopped_reason()
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

bool RaystarNode::buildPathMsg(const PathSolution& solution,
                               const GridMap& grid_map,
                               const std::string& frame_id,
                               size_t max_path_points,
                               nav_msgs::msg::Path& msg,
                               std::string& error) const {
  msg = nav_msgs::msg::Path{};
  error.clear();
  msg.header.stamp = now();
  msg.header.frame_id = frame_id;

  std::vector<Point2d> interpolated;
  if (!interpolateProjectedPath(solution, max_path_points, interpolated, error)) {
    return false;
  }
  msg.poses.reserve(interpolated.size());
  for (const auto& point : interpolated) {
    geometry_msgs::msg::PoseStamped pose;
    pose.header = msg.header;
    pose.pose.orientation.w = 1.0;
    double wx, wy;
    if (!continuousGridToWorld(grid_map, point, wx, wy)) {
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

bool RaystarNode::buildTopologyPathMsg(const PathSolution& solution,
                                       const GridMap& grid_map,
                                       const std::string& frame_id,
                                       nav_msgs::msg::Path& msg,
                                       std::string& error) const {
  msg = nav_msgs::msg::Path{};
  error.clear();
  msg.header.stamp = now();
  msg.header.frame_id = frame_id;

  const auto projected = solution.projectedPath();
  if (projected.size() < 2) {
    error = "topology path must contain both endpoints";
    return false;
  }
  msg.poses.reserve(projected.size());
  for (const auto& point : projected) {
    geometry_msgs::msg::PoseStamped pose;
    pose.header = msg.header;
    pose.pose.orientation.w = 1.0;
    if (!continuousGridToWorld(
          grid_map, point, pose.pose.position.x, pose.pose.position.y)) {
      error = "topology path contains a point outside the finite world transform";
      msg.poses.clear();
      return false;
    }
    msg.poses.emplace_back(std::move(pose));
  }
  return true;
}

void RaystarNode::publishPolyObstacles(const Polymap& polymap,
                                       const GridMap& grid_map,
                                       const std::string& frame_id,
                                       size_t max_marker_bytes) {
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
  for (const auto& ob : polymap.obstacles()) {
    if (emitted_markers >= marker_limit) {
      truncated = true;
      break;
    }
    marker.id++;
    marker.points.clear();
    marker.colors.clear();

    for (auto it = ob.ordered_vertices_.begin(); it != ob.ordered_vertices_.end(); ++it) {
      auto nxt = std::next(it);
      if (nxt == ob.ordered_vertices_.end())
        nxt = ob.ordered_vertices_.begin();

      if (!canAppendCount(emitted_points, 2, point_limit)) {
        truncated = true;
        break;
      }

      double wx1, wy1, wx2, wy2;
      mapToWorld(grid_map,
                 static_cast<unsigned int>(it->first),
                 static_cast<unsigned int>(it->second),
                 wx1,
                 wy1);
      mapToWorld(grid_map,
                 static_cast<unsigned int>(nxt->first),
                 static_cast<unsigned int>(nxt->second),
                 wx2,
                 wy2);

      geometry_msgs::msg::Point p1, p2;
      p1.x = wx1;
      p1.y = wy1;
      p2.x = wx2;
      p2.y = wy2;
      marker.points.push_back(p1);
      marker.points.push_back(p2);

      std_msgs::msg::ColorRGBA c;
      c.r = 1.0;
      c.a = 1.0;
      marker.colors.push_back(c);
      marker.colors.push_back(c);
      emitted_points += 2;
    }
    if (!marker.points.empty()) {
      array.markers.push_back(marker);
      ++emitted_markers;
    }
    if (truncated)
      break;
  }
  (void)truncated;
  poly_obstacle_pub_->publish(array);
}

void RaystarNode::publishCDT(const Polymap& polymap,
                             const GridMap& grid_map,
                             const std::string& frame_id,
                             size_t max_marker_bytes) {
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
  for (const auto& e : cdt_edges) {
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

    if (e.is_constrained) {
      con_marker.points.push_back(p1);
      con_marker.points.push_back(p2);
    } else {
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

void RaystarNode::publishNonHomotopicPaths(const std::vector<PathSolution>& solutions,
                                           const GridMap& grid_map,
                                           const std::string& frame_id,
                                           size_t max_path_points,
                                           size_t max_marker_bytes) {
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

  for (size_t si = 0; si < solutions.size(); ++si) {
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
    if (!interpolateProjectedPath(solutions[si],
                                  std::min(max_path_points, remaining_points),
                                  interpolated,
                                  interpolation_error)) {
      RCLCPP_WARN(get_logger(),
                  "Skipping path visualization because output limits rejected it: %s",
                  interpolation_error.c_str());
      continue;
    }
    for (const auto& point : interpolated) {
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
    if (!marker.points.empty()) {
      array.markers.push_back(marker);
      ++emitted_markers;
    }
  }
  auto snapshot = std::make_shared<MarkerArray>(std::move(array));
  non_homotopic_pub_->publish(*snapshot);
  if (path_visualization_timer_)
    cached_path_visualization_ = std::move(snapshot);
}

void RaystarNode::clearVisualizationsLocked() noexcept {
  std::string frame_id;
  try {
    // Accepted request frames are bounded by kMaxFrameIdBytes.  If even this
    // small copy cannot be made, an empty frame is still safe for DELETEALL.
    frame_id = last_frame_id_;
  } catch (...) {
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
  try {
    const auto clear_snapshot = makeMarkerSnapshot(frame_id, now());
    for (const auto& publisher :
         {non_homotopic_pub_, poly_obstacle_pub_, debug_tree_pub_, cdt_pub_}) {
      try {
        publisher->publish(clear_snapshot);
      } catch (...) {
        // Continue clearing the other durable topics.  The local caches and
        // Core state have already been released.
      }
    }
  } catch (...) {
    // Even the small DELETEALL snapshot could not be allocated.
  }
}

void RaystarNode::republishCachedPathVisualization() {
  std::unique_lock<std::mutex> planner_lock(planner_cache_mutex_, std::try_to_lock);
  if (!planner_lock.owns_lock() || planning_busy_.load(std::memory_order_acquire)) {
    return;
  }
  if (!cached_path_visualization_)
    return;

  // Keep publish inside planner_cache_mutex_.  Copying the shared_ptr and
  // publishing after unlock would allow a new request to publish DELETEALL
  // first and this old snapshot second, resurrecting stale paths.
  try {
    if (non_homotopic_pub_->get_subscription_count() == 0)
      return;
    non_homotopic_pub_->publish(*cached_path_visualization_);
  } catch (const std::exception& exception) {
    RCLCPP_WARN_THROTTLE(get_logger(),
                         *get_clock(),
                         5000,
                         "Could not republish cached path visualization: %s",
                         exception.what());
  } catch (...) {
    RCLCPP_WARN_THROTTLE(
      get_logger(), *get_clock(), 5000, "Could not republish cached path visualization");
  }
}

void RaystarNode::publishDebugTree(const std::vector<raystar::Node>& nodes,
                                   const GridMap& grid_map,
                                   const std::string& frame_id,
                                   size_t max_debug_nodes,
                                   size_t max_marker_bytes) {
  const auto stamp = now();
  auto array = makeMarkerSnapshot(frame_id, stamp);
  const size_t node_count = std::min(nodes.size(), max_debug_nodes);
  if (node_count == 0) {
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
  for (size_t i = 0; i < node_count; ++i) {
    const auto& n = nodes[i];
    const double f = (n.gCost() + n.hCost()) * resolution;
    if (f < min_f)
      min_f = f;
    if (f > max_f)
      max_f = f;
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

  for (size_t i = 0; i < node_count; ++i) {
    if (emitted_markers >= node_marker_limit || emitted_points >= point_limit)
      break;
    const auto& n = nodes[i];
    double wx, wy;
    mapToWorld(grid_map,
               static_cast<unsigned int>(n.seed().first),
               static_cast<unsigned int>(n.seed().second),
               wx,
               wy);

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

    if (n.parentIndex() >= 0 && static_cast<size_t>(n.parentIndex()) < node_count) {
      const auto& parent = nodes[static_cast<size_t>(n.parentIndex())];
      double pwx, pwy;
      mapToWorld(grid_map,
                 static_cast<unsigned int>(parent.seed().first),
                 static_cast<unsigned int>(parent.seed().second),
                 pwx,
                 pwy);
      geometry_msgs::msg::Point pp;
      pp.x = pwx;
      pp.y = pwy;
      if (canAppendCount(emitted_points, 2, point_limit)) {
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
    oss << "N" << n.index() << " G=" << std::fixed << std::setprecision(1) << n.gCost() * resolution
        << " F=" << std::fixed << std::setprecision(1) << f;
    text_marker.text = oss.str();
    if (emitted_markers >= node_marker_limit)
      break;
    array.markers.push_back(text_marker);
    ++emitted_markers;

    if (!n.visibility().empty() && emitted_markers < node_marker_limit &&
        emitted_points < point_limit) {
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
      for (const auto& v : n.visibility()) {
        if (!canAppendCount(emitted_points, 1, point_limit))
          break;
        geometry_msgs::msg::Point vp;
        vp.x = grid_map.origin_x + v.first * grid_map.resolution;
        vp.y = grid_map.origin_y + v.second * grid_map.resolution;
        visreg_marker.points.push_back(vp);
        ++emitted_points;
      }
      if (n.index() > 0) {
        if (canAppendCount(emitted_points, 1, point_limit)) {
          geometry_msgs::msg::Point sp;
          sp.x = grid_map.origin_x + static_cast<double>(n.seed().first) * grid_map.resolution;
          sp.y = grid_map.origin_y + static_cast<double>(n.seed().second) * grid_map.resolution;
          visreg_marker.points.push_back(sp);
          ++emitted_points;
        }
      }
      if (visreg_marker.points.size() > 2 && canAppendCount(emitted_points, 1, point_limit)) {
        visreg_marker.points.push_back(visreg_marker.points.front());
        ++emitted_points;
      }
      if (!visreg_marker.points.empty() && emitted_markers < node_marker_limit) {
        array.markers.push_back(visreg_marker);
        ++emitted_markers;
      }

      if (!n.fullVisibility().empty() && n.fullVisibility().size() != n.visibility().size() &&
          emitted_markers < node_marker_limit && emitted_points < point_limit) {
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
        for (const auto& v : n.fullVisibility()) {
          if (!canAppendCount(emitted_points, 1, point_limit))
            break;
          geometry_msgs::msg::Point vp;
          vp.x = grid_map.origin_x + v.first * grid_map.resolution;
          vp.y = grid_map.origin_y + v.second * grid_map.resolution;
          full_visreg_marker.points.push_back(vp);
          ++emitted_points;
        }
        if (n.index() > 0) {
          if (canAppendCount(emitted_points, 1, point_limit)) {
            geometry_msgs::msg::Point sp;
            sp.x = grid_map.origin_x + static_cast<double>(n.seed().first) * grid_map.resolution;
            sp.y = grid_map.origin_y + static_cast<double>(n.seed().second) * grid_map.resolution;
            full_visreg_marker.points.push_back(sp);
            ++emitted_points;
          }
        }
        if (full_visreg_marker.points.size() > 2 &&
            canAppendCount(emitted_points, 1, point_limit)) {
          full_visreg_marker.points.push_back(full_visreg_marker.points.front());
          ++emitted_points;
        }
        if (!full_visreg_marker.points.empty() && emitted_markers < node_marker_limit) {
          array.markers.push_back(full_visreg_marker);
          ++emitted_markers;
        }
      }

      if (n.localShortestPath().size() >= 2 && emitted_markers < node_marker_limit &&
          emitted_points < point_limit) {
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
        for (const auto& wp : n.localShortestPath()) {
          if (!canAppendCount(emitted_points, 1, point_limit))
            break;
          double tpx, tpy;
          mapToWorld(grid_map,
                     static_cast<unsigned int>(wp.first),
                     static_cast<unsigned int>(wp.second),
                     tpx,
                     tpy);
          geometry_msgs::msg::Point tp;
          tp.x = tpx;
          tp.y = tpy;
          tpath_marker.points.push_back(tp);
          ++emitted_points;
        }
        if (!tpath_marker.points.empty() && emitted_markers < node_marker_limit) {
          array.markers.push_back(tpath_marker);
          ++emitted_markers;
        }
      }
    }
  }

  if (emitted_markers < marker_limit) {
    array.markers.push_back(edge_marker);
    ++emitted_markers;
  }
  if (emitted_markers < marker_limit)
    array.markers.push_back(seed_marker);

  debug_tree_pub_->publish(array);
}

}  // namespace raystar
