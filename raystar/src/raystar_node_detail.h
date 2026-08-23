#ifndef RAYSTAR_NODE_DETAIL_H_
#define RAYSTAR_NODE_DETAIL_H_

// Package-private helpers shared by the RaystarNode translation units
// (raystar_node*.cpp).  Formerly an anonymous namespace inside the single
// raystar_node.cpp.

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

namespace node_impl {


// This tolerance validates quaternion encoding only. Rotation components are
// still required to be exactly zero, so no geometric yaw/tilt is ignored.
constexpr double kQuaternionNormTolerance = 1e-12;
constexpr int64_t kDefaultOccupiedThreshold = 99;
constexpr int64_t kDefaultMaxK = 100;
constexpr int64_t kDefaultMaxCostBoundedPaths =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxCostBoundedPaths);
constexpr int64_t kDefaultMaxMultiGoalCount =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxMultiGoalCount);
constexpr int64_t kDefaultMaxTransitionReferences =
  static_cast<int64_t>(PlanningLimits::kDefaultMaxTransitionReferences);
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
  {"max_transition_references",
   kDefaultMaxTransitionReferences,
   1,
   kMaxIntParameterValue,
   "Maximum flattened rooted references accepted by one UPS batch.",
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

inline const IntegerParameterSpec* findIntegerParameterSpec(const std::string& name) {
  const auto found =
    std::find_if(kIntegerParameterSpecs.begin(),
                 kIntegerParameterSpecs.end(),
                 [&name](const IntegerParameterSpec& spec) { return name == spec.name; });
  return found == kIntegerParameterSpecs.end() ? nullptr : &*found;
}

inline rcl_interfaces::msg::ParameterDescriptor makeParameterDescriptor(const IntegerParameterSpec& spec) {
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

inline rcl_interfaces::msg::SetParametersResult validateParameterChanges(
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

inline size_t markerTopicBudget(size_t total_budget) {
  return std::max<size_t>(1, total_budget / kVisualizationTopicCount);
}

inline size_t markerPointLimit(size_t topic_budget) {
  return topic_budget / kMarkerPointBytes;
}

inline size_t markerEntryLimit(size_t topic_budget) {
  return topic_budget / kMarkerEntryBytes;
}

inline bool canAppendCount(size_t current, size_t addition, size_t limit) {
  return addition <= limit && current <= limit - addition;
}

inline visualization_msgs::msg::MarkerArray makeMarkerSnapshot(const std::string& frame_id,
                                                        const rclcpp::Time& stamp) {
  visualization_msgs::msg::MarkerArray array;
  visualization_msgs::msg::Marker clear_marker;
  clear_marker.header.frame_id = frame_id;
  clear_marker.header.stamp = stamp;
  clear_marker.action = visualization_msgs::msg::Marker::DELETEALL;
  array.markers.emplace_back(std::move(clear_marker));
  return array;
}

inline bool validatePlanarPose(const geometry_msgs::msg::PoseStamped& pose,
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
inline bool continuousGridToWorld(const GridMap& grid_map, const Point2d& point, double& wx, double& wy) {
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
inline bool countInterpolatedPathPoints(const PathSolution& solution,
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

inline bool interpolateProjectedPath(const PathSolution& solution,
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

inline bool checkedAdd(size_t& total, size_t addition, size_t limit) {
  if (addition > limit || total > limit - addition)
    return false;
  total += addition;
  return true;
}

inline bool estimatePathResponseBytes(size_t point_count, size_t frame_id_bytes, size_t& bytes) {
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

inline bool publishedPathLength(const nav_msgs::msg::Path& path,
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
    if (!std::isfinite(segment) || !certificate.addSegment(first.x, first.y, second.x, second.y)) {
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

inline bool publishedLengthsMatch(double first, double second) {
  const double tolerance =
    std::max(kPublishedLengthAbsoluteToleranceMeters,
             kPublishedLengthRelativeTolerance * std::max(std::abs(first), std::abs(second)));
  return std::abs(first - second) <= tolerance;
}

enum class MetricBoundEligibility { within_bound, outside_bound, invalid };

inline MetricBoundEligibility classifyPathViewMetricBound(const BoundedPathView& path,
                                                   const GridMap& map,
                                                   double inclusive_bound,
                                                   const StopToken& stop_token) {
  if (!std::isfinite(path.path_cost) || path.path_cost < 0.0 || !std::isfinite(inclusive_bound) ||
      inclusive_bound < 0.0) {
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
    if (have_previous && !certificate.addSegment(previous_x, previous_y, world_x, world_y)) {
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
          {static_cast<double>(turning_point.first), static_cast<double>(turning_point.second)})) {
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
  if (!std::isfinite(core_metric_cost) || !publishedLengthsMatch(core_metric_cost, upper_bound)) {
    return MetricBoundEligibility::invalid;
  }
  return upper_bound <= inclusive_bound ? MetricBoundEligibility::within_bound
                                        : MetricBoundEligibility::outside_bound;
}

inline MetricBoundEligibility classifyTopologyMetricBound(const nav_msgs::msg::Path& topology_path,
                                                   double core_cost,
                                                   double inclusive_bound,
                                                   const StopPredicate& stop_requested = {}) {
  if (!std::isfinite(core_cost) || core_cost < 0.0 || !std::isfinite(inclusive_bound) ||
      inclusive_bound < 0.0) {
    return MetricBoundEligibility::invalid;
  }
  double rounded_length = std::numeric_limits<double>::quiet_NaN();
  double upper_bound = std::numeric_limits<double>::quiet_NaN();
  if (!publishedPathLength(topology_path, rounded_length, upper_bound, stop_requested) ||
      !publishedLengthsMatch(core_cost, upper_bound)) {
    return MetricBoundEligibility::invalid;
  }
  return upper_bound <= inclusive_bound ? MetricBoundEligibility::within_bound
                                        : MetricBoundEligibility::outside_bound;
}

inline bool finalizePublishedPathResult(raystar_interfaces::msg::PathResult& result,
                                 double core_cost,
                                 double inclusive_bound = std::numeric_limits<double>::infinity(),
                                 const StopPredicate& stop_requested = {}) {
  if (!std::isfinite(core_cost) || core_cost < 0.0 || std::isnan(inclusive_bound))
    return false;
  double dense_length = std::numeric_limits<double>::quiet_NaN();
  double topology_length = std::numeric_limits<double>::quiet_NaN();
  double dense_upper_bound = std::numeric_limits<double>::quiet_NaN();
  double topology_upper_bound = std::numeric_limits<double>::quiet_NaN();
  if (!publishedPathLength(result.path, dense_length, dense_upper_bound, stop_requested) ||
      !publishedPathLength(
        result.topology_path, topology_length, topology_upper_bound, stop_requested)) {
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

inline bool finalizeCostBoundedPublishedPathResult(raystar_interfaces::msg::PathResult& result,
                                            double core_cost,
                                            double inclusive_bound,
                                            const StopPredicate& stop_requested = {}) {
  if (finalizePublishedPathResult(result, core_cost, inclusive_bound, stop_requested))
    return true;

  // Interpolation is visualization-only. Rounding its intermediate world
  // coordinates can make the serialized dense polyline microscopically
  // longer than the exact topology polyline. At an inclusive one-ULP bound,
  // retry with the collision-free topology poses themselves instead of
  // discarding a mathematically eligible path or weakening the certificate.
  result.path = result.topology_path;
  return finalizePublishedPathResult(result, core_cost, inclusive_bound, stop_requested);
}

inline bool hasMetricFaithfulWorldTransform(
  double origin_x, double origin_y, double extent_x, double extent_y, double resolution) {
  // The published path contract allows a small representation discrepancy
  // from Core's metric cost.  Admit only transforms whose binary64 quantum is
  // comfortably below that budget at both the multiply and add magnitudes.
  // Per-path validation below remains the authoritative final check.
  const double representation_budget = std::max(kPublishedLengthAbsoluteToleranceMeters,
                                                kPublishedLengthRelativeTolerance * resolution) /
                                       32.0;
  // The representation tolerance alone is insufficient for very small map
  // resolutions: a large origin can have an ULP smaller than 1e-9 metres yet
  // still collapse many adjacent cells. Requiring the coordinate quantum to
  // be much smaller than one cell makes every rounded adjacent transform
  // strictly monotone as well as metrically distinguishable.
  const double coordinate_budget = std::min(representation_budget, resolution / 32.0);
  const double max_world_x = std::fma(1.0, extent_x, origin_x);
  const double max_world_y = std::fma(1.0, extent_y, origin_y);
  const double first_world_x = std::fma(1.0, resolution, origin_x);
  const double first_world_y = std::fma(1.0, resolution, origin_y);
  const std::array<double, 8> transform_values{
    origin_x, origin_y, extent_x, extent_y, first_world_x, first_world_y, max_world_x, max_world_y};
  if (!std::isfinite(coordinate_budget) || coordinate_budget <= 0.0 ||
      !std::all_of(transform_values.begin(), transform_values.end(), [&](double value) {
        return std::isfinite(value) && binary64SpacingForSearchPadding(value) <= coordinate_budget;
      })) {
    return false;
  }
  if (!(first_world_x > origin_x) || !(first_world_y > origin_y) || !(max_world_x > origin_x) ||
      !(max_world_y > origin_y)) {
    return false;
  }
  const auto distance_is_represented = [](double actual, double expected) {
    if (!std::isfinite(actual) || actual <= 0.0)
      return false;
    const long double error =
      std::abs(static_cast<long double>(actual) - static_cast<long double>(expected));
    const long double tolerance =
      std::max(static_cast<long double>(kPublishedLengthAbsoluteToleranceMeters),
               static_cast<long double>(kPublishedLengthRelativeTolerance) *
                 static_cast<long double>(expected));
    return error <= tolerance;
  };
  return distance_is_represented(first_world_x - origin_x, resolution) &&
         distance_is_represented(first_world_y - origin_y, resolution) &&
         distance_is_represented(max_world_x - origin_x, extent_x) &&
         distance_is_represented(max_world_y - origin_y, extent_y);
}

inline bool appendResponseNotice(std::string& message,
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

inline void setBoundedExceptionMessage(std::string& message, const char* detail) {
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

inline uint32_t boundedUint32(size_t value) {
  return static_cast<uint32_t>(
    std::min(value, static_cast<size_t>(std::numeric_limits<uint32_t>::max())));
}

inline uint64_t boundedUint64(size_t value) {
  if constexpr (sizeof(size_t) > sizeof(uint64_t)) {
    return static_cast<uint64_t>(
      std::min(value, static_cast<size_t>(std::numeric_limits<uint64_t>::max())));
  }
  return static_cast<uint64_t>(value);
}

inline uint16_t planningLimitMask(PlanningLimitReached limit) {
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

inline void markPathOutputIncomplete(PlanningResultInfo& info,
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


}  // namespace node_impl

}  // namespace raystar

#endif  // RAYSTAR_NODE_DETAIL_H_
