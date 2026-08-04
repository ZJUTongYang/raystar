#pragma once

#include <raystar/coordinate_utils.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

namespace raystar {

// Internal, header-defined helpers keep the metric-to-Core search-superset
// proof directly unit-testable without exposing it as public package API.
inline double binary64SpacingForSearchPadding(double value) {
  if (!std::isfinite(value))
    return std::numeric_limits<double>::infinity();
  const double upward = std::nextafter(value, std::numeric_limits<double>::infinity()) - value;
  const double downward = value - std::nextafter(value, -std::numeric_limits<double>::infinity());
  return std::max(upward, downward);
}

inline double saturatedUpwardPositiveProduct(double first, double second) {
  const double rounded = first * second;
  if (!std::isfinite(rounded))
    return std::numeric_limits<double>::max();
  const double successor = std::nextafter(rounded, std::numeric_limits<double>::infinity());
  // FE_DOWNWARD and FE_TOWARDZERO can round a positive overflowing product
  // to DBL_MAX. There is no finite successor in that case, but DBL_MAX is
  // still the correct search-superset saturation because every Core path
  // cost is finite.
  return std::isfinite(successor) ? successor : std::numeric_limits<double>::max();
}

inline double saturatedUpwardPositiveSum(double first, double second) {
  const double rounded = first + second;
  if (!std::isfinite(rounded))
    return std::numeric_limits<double>::max();
  const double successor = std::nextafter(rounded, std::numeric_limits<double>::infinity());
  // As above, saturate instead of turning a finite request into an invalid
  // infinite padded bound under a directed process rounding mode.
  return std::isfinite(successor) ? successor : std::numeric_limits<double>::max();
}

inline bool paddedMetricBoundForGridSearch(const GridMap& map,
                                           double metric_bound,
                                           size_t maximum_search_nodes,
                                           double& padded_bound) {
  padded_bound = std::numeric_limits<double>::quiet_NaN();
  if (!hasValidWorldTransform(map) || !std::isfinite(metric_bound) || metric_bound <= 0.0 ||
      maximum_search_nodes == 0) {
    return false;
  }

  const double resolution = static_cast<double>(map.resolution);
  const double max_world_x = std::fma(static_cast<double>(map.width), resolution, map.origin_x);
  const double max_world_y = std::fma(static_cast<double>(map.height), resolution, map.origin_y);
  const std::array<double, 4> operation_values{
    map.origin_x, map.origin_y, max_world_x, max_world_y};
  double operation_spacing = 0.0;
  for (const double value : operation_values) {
    if (!std::isfinite(value))
      return false;
    operation_spacing = std::max(operation_spacing, binary64SpacingForSearchPadding(value));
  }
  if (!std::isfinite(operation_spacing))
    return false;

  // continuousGridToWorld uses one explicitly rounded fma per coordinate.
  // It differs from the exact affine transform by at most one local spacing
  // under every IEEE rounding mode. The reverse triangle inequality bounds
  // one segment's possible serialized-world contraction by
  // 2 endpoints * sqrt(2) * spacing.
  constexpr double kSquareRootTwoUpper = 0x1.6a09e667f3bcdp+0;
  double segment_contraction = saturatedUpwardPositiveProduct(operation_spacing, 2.0);
  segment_contraction = saturatedUpwardPositiveProduct(segment_contraction, kSquareRootTwoUpper);

  // The root node can connect a one-segment path; every later search node can
  // add one turning point and hence one segment. When max_nodes=N nodes have
  // already been constructed, the queue may still contain a child optimistic
  // certificate with N turns plus its final goal segment. Its N+1 segments
  // must remain in the padded search superset so the next loop iteration
  // reports LIMIT_MAX_NODES instead of incorrectly claiming bound exhaustion.
  // max_path_points is only an output/retention budget and cannot bound such
  // an undiscovered frontier path.
  const size_t maximum_segment_count = maximum_search_nodes == std::numeric_limits<size_t>::max()
                                         ? maximum_search_nodes
                                         : maximum_search_nodes + 1u;
  double segment_count_upper = static_cast<double>(maximum_segment_count);
  if (!std::isfinite(segment_count_upper)) {
    segment_count_upper = std::numeric_limits<double>::max();
  } else {
    const double successor =
      std::nextafter(segment_count_upper, std::numeric_limits<double>::infinity());
    segment_count_upper = std::isfinite(successor) ? successor : std::numeric_limits<double>::max();
  }
  const double path_contraction =
    saturatedUpwardPositiveProduct(segment_contraction, segment_count_upper);
  padded_bound = saturatedUpwardPositiveSum(metric_bound, path_contraction);
  return std::isfinite(padded_bound) && padded_bound > 0.0;
}

}  // namespace raystar
