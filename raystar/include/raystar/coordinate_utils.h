#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include <raystar/cooperative_stop.h>

namespace raystar {

// C++17 has no std::numbers::pi.  Keep one portable project constant instead
// of relying on non-standard libc pi macros.
inline constexpr double kPi = 3.141592653589793238462643383279502884;
inline constexpr double kTwoPi = 2.0 * kPi;

/// Grid map data structure (binary occupancy: 0=free, 1=occupied)
struct GridMap {
  std::vector<uint8_t> data;
  unsigned int width = 0;
  unsigned int height = 0;
  float resolution = 0.05f;
  double origin_x = 0.0;
  double origin_y = 0.0;

  /// Access cell value (row-major). Returns 1 (occupied) for out-of-bounds.
  inline uint8_t operator()(unsigned int x, unsigned int y) const {
    if (x >= width || y >= height)
      return 1;
    return data[y * width + x];
  }

  inline uint8_t& operator()(unsigned int x, unsigned int y) {
    return data[y * width + x];
  }

  inline uint8_t at(unsigned int x, unsigned int y) const {
    if (x >= width || y >= height)
      return 1;
    return data[y * width + x];
  }
};

/// A point in the map's continuous grid coordinate system together with the
/// integer cell selected by floor().  The cell is deliberately kept separate
/// from the point: a cell is used for occupancy/indexing while the point is
/// used by the geometric planner.
struct ContinuousGridPoint {
  double x = 0.0;
  double y = 0.0;
  unsigned int cell_x = 0;
  unsigned int cell_y = 0;
};

/// Return true when the unrotated world/grid transform can be evaluated
/// without producing non-finite coordinates.
inline bool hasValidWorldTransform(const GridMap& map) {
  return std::isfinite(static_cast<double>(map.resolution)) && map.resolution > 0.0f &&
         std::isfinite(map.origin_x) && std::isfinite(map.origin_y);
}

// World/grid conversion can lose a couple of low bits when a caller first
// computes `origin + integer * resolution` and we subsequently subtract the
// origin and divide by the (float-stored) resolution.  In that case a point
// which is exactly a grid boundary may arrive just below the integer and be
// assigned to the preceding cell by floor().  Recover only boundaries that
// round-trip to the supplied world value (within two ULPs); this is narrow
// enough to preserve genuinely sub-cell positions while making the two
// conversion directions agree at representable boundaries.
inline double canonicalizeWorldGridCoordinate(
  double world, double origin, double resolution, unsigned int extent, double converted) {
  if (!std::isfinite(converted))
    return converted;

  const double nearest = std::round(converted);
  if (!std::isfinite(nearest) || nearest < 0.0 || nearest > static_cast<double>(extent)) {
    return converted;
  }

  const double reconstructed = origin + nearest * resolution;
  if (!std::isfinite(reconstructed))
    return converted;

  const double upward_ulp = std::abs(
    std::nextafter(reconstructed, std::numeric_limits<double>::infinity()) - reconstructed);
  const double downward_ulp = std::abs(
    reconstructed - std::nextafter(reconstructed, -std::numeric_limits<double>::infinity()));
  const double tolerance = 2.0 * std::max(upward_ulp, downward_ulp);
  if (world == reconstructed ||
      (std::isfinite(tolerance) && std::abs(world - reconstructed) <= tolerance)) {
    return nearest;
  }
  return converted;
}

/// World coordinates (meters) -> Grid coordinates (cells).
/// Returns false if outside the map.
inline bool worldToMap(
  const GridMap& map, double wx, double wy, unsigned int& mx, unsigned int& my) {
  if (!hasValidWorldTransform(map) || map.width == 0 || map.height == 0 || !std::isfinite(wx) ||
      !std::isfinite(wy)) {
    return false;
  }

  double dx = (wx - map.origin_x) / map.resolution;
  double dy = (wy - map.origin_y) / map.resolution;
  dx = canonicalizeWorldGridCoordinate(
    wx, map.origin_x, static_cast<double>(map.resolution), map.width, dx);
  dy = canonicalizeWorldGridCoordinate(
    wy, map.origin_y, static_cast<double>(map.resolution), map.height, dy);
  if (!std::isfinite(dx) || !std::isfinite(dy) || dx < 0.0 || dy < 0.0 ||
      dx >= static_cast<double>(map.width) || dy >= static_cast<double>(map.height)) {
    return false;
  }

  const unsigned int converted_x = static_cast<unsigned int>(dx);
  const unsigned int converted_y = static_cast<unsigned int>(dy);
  mx = converted_x;
  my = converted_y;
  return true;
}

/// Convert world coordinates to continuous grid coordinates without losing
/// the sub-cell position.  `cell_x/cell_y` are exactly floor(x/y), and all
/// output arguments are left unchanged when conversion fails.
inline bool worldToContinuousMap(const GridMap& map,
                                 double wx,
                                 double wy,
                                 double& gx,
                                 double& gy,
                                 unsigned int& cell_x,
                                 unsigned int& cell_y) {
  if (!hasValidWorldTransform(map) || map.width == 0 || map.height == 0 || !std::isfinite(wx) ||
      !std::isfinite(wy)) {
    return false;
  }

  double converted_x = (wx - map.origin_x) / static_cast<double>(map.resolution);
  double converted_y = (wy - map.origin_y) / static_cast<double>(map.resolution);
  converted_x = canonicalizeWorldGridCoordinate(
    wx, map.origin_x, static_cast<double>(map.resolution), map.width, converted_x);
  converted_y = canonicalizeWorldGridCoordinate(
    wy, map.origin_y, static_cast<double>(map.resolution), map.height, converted_y);
  if (!std::isfinite(converted_x) || !std::isfinite(converted_y) || converted_x < 0.0 ||
      converted_y < 0.0 || converted_x >= static_cast<double>(map.width) ||
      converted_y >= static_cast<double>(map.height)) {
    return false;
  }

  const double floored_x = std::floor(converted_x);
  const double floored_y = std::floor(converted_y);
  if (floored_x < 0.0 || floored_y < 0.0 || floored_x >= static_cast<double>(map.width) ||
      floored_y >= static_cast<double>(map.height) ||
      floored_x > static_cast<double>(std::numeric_limits<unsigned int>::max()) ||
      floored_y > static_cast<double>(std::numeric_limits<unsigned int>::max())) {
    return false;
  }

  gx = converted_x;
  gy = converted_y;
  cell_x = static_cast<unsigned int>(floored_x);
  cell_y = static_cast<unsigned int>(floored_y);
  return true;
}

/// Convenience overload returning the continuous point and its cell in one
/// value.  It is useful at API boundaries where keeping the two values paired
/// prevents accidentally replacing the point with its cell.
inline bool worldToContinuousMap(const GridMap& map,
                                 double wx,
                                 double wy,
                                 ContinuousGridPoint& result) {
  ContinuousGridPoint candidate;
  if (!worldToContinuousMap(
        map, wx, wy, candidate.x, candidate.y, candidate.cell_x, candidate.cell_y)) {
    return false;
  }
  result = candidate;
  return true;
}

/// Grid coordinates (cells) -> World coordinates (meters).
inline bool mapToWorld(
  const GridMap& map, unsigned int mx, unsigned int my, double& wx, double& wy) {
  wx = std::numeric_limits<double>::quiet_NaN();
  wy = std::numeric_limits<double>::quiet_NaN();
  if (!hasValidWorldTransform(map))
    return false;

  wx = map.origin_x + static_cast<double>(mx) * map.resolution;
  wy = map.origin_y + static_cast<double>(my) * map.resolution;
  return std::isfinite(wx) && std::isfinite(wy);
}

/// Convert a continuous grid point to world coordinates.  Unlike the legacy
/// integer mapToWorld() overload this function never rounds or indexes a cell.
inline bool continuousMapToWorld(const GridMap& map, double gx, double gy, double& wx, double& wy) {
  wx = std::numeric_limits<double>::quiet_NaN();
  wy = std::numeric_limits<double>::quiet_NaN();
  if (!hasValidWorldTransform(map) || !std::isfinite(gx) || !std::isfinite(gy) || gx < 0.0 ||
      gy < 0.0 || gx >= static_cast<double>(map.width) || gy >= static_cast<double>(map.height)) {
    return false;
  }

  const double converted_x = map.origin_x + gx * static_cast<double>(map.resolution);
  const double converted_y = map.origin_y + gy * static_cast<double>(map.resolution);
  if (!std::isfinite(converted_x) || !std::isfinite(converted_y))
    return false;
  wx = converted_x;
  wy = converted_y;
  return true;
}

// Alternate spellings retained as small, header-only aliases for callers that
// use the conventional "continuous" suffix.
inline bool worldToMapContinuous(const GridMap& map,
                                 double wx,
                                 double wy,
                                 double& gx,
                                 double& gy,
                                 unsigned int& cell_x,
                                 unsigned int& cell_y) {
  return worldToContinuousMap(map, wx, wy, gx, gy, cell_x, cell_y);
}

inline bool mapContinuousToWorld(const GridMap& map, double gx, double gy, double& wx, double& wy) {
  return continuousMapToWorld(map, gx, gy, wx, wy);
}

/// Normalize a finite angle to (-pi, pi]. Non-finite input returns NaN.
inline double normalize_angle(double angle) {
  if (!std::isfinite(angle))
    return std::numeric_limits<double>::quiet_NaN();

  double normalized = std::remainder(angle, kTwoPi);
  if (normalized <= -kPi)
    normalized += kTwoPi;
  // Canonicalize negative zero as well as exact multiples of two pi.
  return normalized == 0.0 ? 0.0 : normalized;
}

/// Normalize a finite angle to [0, 2*pi). Non-finite input returns NaN.
inline double normalize_angle_positive(double angle) {
  if (!std::isfinite(angle))
    return std::numeric_limits<double>::quiet_NaN();

  double normalized = std::fmod(angle, kTwoPi);
  if (normalized < 0.0)
    normalized += kTwoPi;
  // Adding two pi to a tiny negative remainder can round to exactly two pi.
  // Zero is the canonical representable value for that wrap point.
  if (normalized >= kTwoPi || normalized == 0.0)
    return 0.0;
  return normalized;
}

/// Return true when `angle` lies in the counter-clockwise sweep from
/// `start_angle` to `end_angle`. The interval may cross -PI/PI and may span
/// more than PI; both comparisons use the same [0, 2*PI) relative-angle
/// representation.
inline bool isAngleWithinCounterClockwiseSweep(double angle,
                                               double start_angle,
                                               double end_angle,
                                               double epsilon = 1e-7) {
  if (!std::isfinite(angle) || !std::isfinite(start_angle) || !std::isfinite(end_angle) ||
      !std::isfinite(epsilon)) {
    return false;
  }

  const double span = normalize_angle_positive(end_angle - start_angle);
  const double relative_angle = normalize_angle_positive(angle - start_angle);
  const double tolerance = std::max(0.0, epsilon);
  return relative_angle <= span + tolerance || kTwoPi - relative_angle <= tolerance;
}

/// Orientation of triplet (p, q, r). Returns 0=collinear, 1=clockwise, 2=counterclockwise.
inline int orientation(std::pair<int, int> p, std::pair<int, int> q, std::pair<int, int> r) {
  // Each int difference needs 33 signed bits, and the determinant can need 65
  // bits.  Compute the two signed products as uint64 magnitudes and compare
  // them without ever materializing their potentially overflowing difference.
  struct SignedProduct {
    std::uint64_t magnitude = 0;
    bool negative = false;
  };

  const auto difference = [](int lhs, int rhs) {
    return static_cast<std::int64_t>(lhs) - static_cast<std::int64_t>(rhs);
  };
  const auto unsigned_magnitude = [](std::int64_t value) {
    return value < 0 ? static_cast<std::uint64_t>(-(value + 1)) + 1u
                     : static_cast<std::uint64_t>(value);
  };
  const auto product = [&](std::int64_t lhs, std::int64_t rhs) {
    SignedProduct result;
    result.magnitude = unsigned_magnitude(lhs) * unsigned_magnitude(rhs);
    result.negative = result.magnitude != 0 && ((lhs < 0) != (rhs < 0));
    return result;
  };

  const SignedProduct first = product(difference(q.second, p.second), difference(r.first, q.first));
  const SignedProduct second =
    product(difference(q.first, p.first), difference(r.second, q.second));

  int determinant_sign = 0;
  if (first.negative != second.negative) {
    determinant_sign = first.negative ? -1 : 1;
  } else if (first.magnitude != second.magnitude) {
    const bool first_has_larger_magnitude = first.magnitude > second.magnitude;
    determinant_sign = first.negative ? (first_has_larger_magnitude ? -1 : 1)
                                      : (first_has_larger_magnitude ? 1 : -1);
  }

  if (determinant_sign == 0)
    return 0;
  return determinant_sign > 0 ? 1 : 2;
}

/// Check if point r lies on segment pq. PRECONDITION: r is collinear with p, q.
inline bool onSegment(std::pair<int, int> p, std::pair<int, int> q, std::pair<int, int> r) {
  return r.first >= std::min(p.first, q.first) && r.first <= std::max(p.first, q.first) &&
         r.second >= std::min(p.second, q.second) && r.second <= std::max(p.second, q.second);
}

inline bool segmentsIntersect(std::pair<int, int> p1,
                              std::pair<int, int> p2,
                              std::pair<int, int> p3,
                              std::pair<int, int> p4) {
  int o1 = orientation(p1, p2, p3);
  int o2 = orientation(p1, p2, p4);
  int o3 = orientation(p3, p4, p1);
  int o4 = orientation(p3, p4, p2);

  if (o1 != o2 && o3 != o4)
    return true;

  if (o1 == 0 && onSegment(p1, p2, p3))
    return true;
  if (o2 == 0 && onSegment(p1, p2, p4))
    return true;
  if (o3 == 0 && onSegment(p3, p4, p1))
    return true;
  if (o4 == 0 && onSegment(p3, p4, p2))
    return true;

  return false;
}

/// Cooperative variant of newPathSelfCrosses(). The result is committed only
/// after the complete check succeeds; on stop, self_crosses is left unchanged.
inline OperationStatus newPathSelfCrosses(const std::vector<std::pair<int, int>>& path,
                                          std::pair<int, int> new_waypoint,
                                          bool& self_crosses,
                                          const StopToken& stop_token) {
  if (stop_token.poll())
    return OperationStatus::stopped;

  if (path.size() < 2) {
    self_crosses = false;
    return OperationStatus::success;
  }

  const auto new_start = path.back();
  const auto new_end = new_waypoint;

  for (size_t i = 0; i + 1 < path.size(); ++i) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    if (path[i] == new_waypoint) {
      if (stop_token.poll())
        return OperationStatus::stopped;
      self_crosses = true;
      return OperationStatus::success;
    }
  }

  for (size_t i = 0; i + 2 < path.size(); ++i) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    if (segmentsIntersect(path[i], path[i + 1], new_start, new_end)) {
      if (stop_token.poll())
        return OperationStatus::stopped;
      self_crosses = true;
      return OperationStatus::success;
    }
  }

  if (stop_token.poll())
    return OperationStatus::stopped;

  self_crosses = false;
  return OperationStatus::success;
}

/// Check if extending a waypoint path with new_waypoint causes self-crossing.
/// Detects vertex revisit (new_waypoint == any existing waypoint except path.back())
/// and geometric segment crossing (new segment intersects any existing segment).
inline bool newPathSelfCrosses(const std::vector<std::pair<int, int>>& path,
                               std::pair<int, int> new_waypoint) {
  bool self_crosses = false;
  const StopToken never_stop;
  (void)newPathSelfCrosses(path, new_waypoint, self_crosses, never_stop);
  return self_crosses;
}

}  // namespace raystar
