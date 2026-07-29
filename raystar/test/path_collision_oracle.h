#ifndef RAYSTAR_TEST_PATH_COLLISION_ORACLE_H_
#define RAYSTAR_TEST_PATH_COLLISION_ORACLE_H_

#include <geometry_msgs/msg/point.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <nav_msgs/msg/path.hpp>

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// This file deliberately implements test-only obstacle-collision and
// self-intersection oracles without including or calling Raystar's coordinate,
// contour, visibility, collision, or pruning helpers.  Keeping the
// implementation independent makes it useful as a regression oracle for final
// paths returned by the ROS API.
namespace raystar::test_oracle
{

enum class ValidationStatus
{
  kCollisionFree,
  kCollision,
  kInconclusive,
};

struct OracleOptions
{
  int occupied_threshold = 99;
  bool allow_unknown = false;
};

struct ValidationResult
{
  static constexpr std::size_t kNoSegment =
    std::numeric_limits<std::size_t>::max();

  ValidationStatus status = ValidationStatus::kInconclusive;
  std::size_t segment_index = kNoSegment;
  long double grid_x = std::numeric_limits<long double>::quiet_NaN();
  long double grid_y = std::numeric_limits<long double>::quiet_NaN();
  std::string diagnostic;

  [[nodiscard]] bool collisionFree() const
  {
    return status == ValidationStatus::kCollisionFree;
  }

  [[nodiscard]] bool collisionDetected() const
  {
    return status == ValidationStatus::kCollision;
  }

  [[nodiscard]] bool conclusive() const
  {
    return status != ValidationStatus::kInconclusive;
  }
};

// Product-level self-intersection semantics for a returned open path.  Two
// consecutive segments are allowed to meet at their construction waypoint;
// every intersection between non-consecutive closed segments is rejected.
// Keeping this result separate from ValidationResult lets tests distinguish
// obstacle collisions from path-topology violations.
enum class SelfIntersectionStatus
{
  kIntersectionFree,
  kSelfIntersection,
  kInconclusive,
};

enum class SelfIntersectionKind
{
  kNone,
  kProperCrossing,
  kNonAdjacentTouch,
  kRepeatedVertex,
  kCollinearOverlap,
};

struct SelfIntersectionResult
{
  static constexpr std::size_t kNoSegment =
    std::numeric_limits<std::size_t>::max();

  SelfIntersectionStatus status = SelfIntersectionStatus::kInconclusive;
  SelfIntersectionKind kind = SelfIntersectionKind::kNone;
  std::size_t first_segment_index = kNoSegment;
  std::size_t second_segment_index = kNoSegment;
  std::string diagnostic;

  [[nodiscard]] bool intersectionFree() const
  {
    return status == SelfIntersectionStatus::kIntersectionFree;
  }

  [[nodiscard]] bool selfIntersectionDetected() const
  {
    return status == SelfIntersectionStatus::kSelfIntersection;
  }

  [[nodiscard]] bool conclusive() const
  {
    return status != SelfIntersectionStatus::kInconclusive;
  }
};

namespace detail
{

// An IEEE-754 double is an exact dyadic rational.  Decode it into an integer
// significand times a power of two, then evaluate orientations with cpp_int.
// This is deliberately independent of the production CGAL predicate and does
// not need an epsilon for proper-crossing, touch, or collinearity decisions.
struct ExactDyadic
{
  boost::multiprecision::cpp_int significand{0};
  int exponent = 0;
};

inline ExactDyadic exactDyadic(double value)
{
  if (value == 0.0)
    return {};

  int exponent = 0;
  const double fraction = std::frexp(value, &exponent);
  constexpr int digits = std::numeric_limits<double>::digits;
  const double scaled = std::ldexp(fraction, digits);
  return ExactDyadic{
    boost::multiprecision::cpp_int(
      static_cast<std::int64_t>(scaled)),
    exponent - digits};
}

inline ExactDyadic subtractExact(
  const ExactDyadic & lhs, const ExactDyadic & rhs)
{
  if (lhs.significand == 0)
    return ExactDyadic{-rhs.significand, rhs.exponent};
  if (rhs.significand == 0)
    return lhs;

  const int exponent = std::min(lhs.exponent, rhs.exponent);
  return ExactDyadic{
    (lhs.significand << (lhs.exponent - exponent)) -
      (rhs.significand << (rhs.exponent - exponent)),
    exponent};
}

inline ExactDyadic multiplyExact(
  const ExactDyadic & lhs, const ExactDyadic & rhs)
{
  if (lhs.significand == 0 || rhs.significand == 0)
    return {};
  return ExactDyadic{
    lhs.significand * rhs.significand, lhs.exponent + rhs.exponent};
}

inline int exactSign(const ExactDyadic & value)
{
  if (value.significand < 0)
    return -1;
  if (value.significand > 0)
    return 1;
  return 0;
}

struct ExactPlanarPoint
{
  double x = 0.0;
  double y = 0.0;
  ExactDyadic exact_x;
  ExactDyadic exact_y;
};

inline ExactPlanarPoint makeExactPlanarPoint(double x, double y)
{
  return ExactPlanarPoint{x, y, exactDyadic(x), exactDyadic(y)};
}

inline bool pointsEqual(
  const ExactPlanarPoint & lhs, const ExactPlanarPoint & rhs)
{
  return lhs.x == rhs.x && lhs.y == rhs.y;
}

inline int exactOrientation(
  const ExactPlanarPoint & first, const ExactPlanarPoint & second,
  const ExactPlanarPoint & third)
{
  const ExactDyadic first_x = subtractExact(second.exact_x, first.exact_x);
  const ExactDyadic first_y = subtractExact(second.exact_y, first.exact_y);
  const ExactDyadic second_x = subtractExact(third.exact_x, first.exact_x);
  const ExactDyadic second_y = subtractExact(third.exact_y, first.exact_y);
  return exactSign(subtractExact(
    multiplyExact(first_x, second_y), multiplyExact(first_y, second_x)));
}

inline bool onClosedSegment(
  const ExactPlanarPoint & first, const ExactPlanarPoint & second,
  const ExactPlanarPoint & point)
{
  return point.x >= std::min(first.x, second.x) &&
    point.x <= std::max(first.x, second.x) &&
    point.y >= std::min(first.y, second.y) &&
    point.y <= std::max(first.y, second.y);
}

inline bool segmentsShareEndpoint(
  const ExactPlanarPoint & first_start,
  const ExactPlanarPoint & first_end,
  const ExactPlanarPoint & second_start,
  const ExactPlanarPoint & second_end)
{
  return pointsEqual(first_start, second_start) ||
    pointsEqual(first_start, second_end) ||
    pointsEqual(first_end, second_start) ||
    pointsEqual(first_end, second_end);
}

inline SelfIntersectionKind classifyClosedSegmentIntersection(
  const ExactPlanarPoint & first_start,
  const ExactPlanarPoint & first_end,
  const ExactPlanarPoint & second_start,
  const ExactPlanarPoint & second_end)
{
  const int first_first = exactOrientation(
    first_start, first_end, second_start);
  const int first_second = exactOrientation(
    first_start, first_end, second_end);
  const int second_first = exactOrientation(
    second_start, second_end, first_start);
  const int second_second = exactOrientation(
    second_start, second_end, first_end);

  if (first_first == 0 && first_second == 0 && second_first == 0 &&
      second_second == 0)
  {
    const bool use_x = first_start.x != first_end.x ||
      second_start.x != second_end.x;
    const double first_lower = use_x ?
      std::min(first_start.x, first_end.x) :
      std::min(first_start.y, first_end.y);
    const double first_upper = use_x ?
      std::max(first_start.x, first_end.x) :
      std::max(first_start.y, first_end.y);
    const double second_lower = use_x ?
      std::min(second_start.x, second_end.x) :
      std::min(second_start.y, second_end.y);
    const double second_upper = use_x ?
      std::max(second_start.x, second_end.x) :
      std::max(second_start.y, second_end.y);
    const double overlap_lower = std::max(first_lower, second_lower);
    const double overlap_upper = std::min(first_upper, second_upper);
    if (overlap_lower > overlap_upper)
      return SelfIntersectionKind::kNone;
    if (overlap_lower < overlap_upper)
      return SelfIntersectionKind::kCollinearOverlap;
    return segmentsShareEndpoint(
      first_start, first_end, second_start, second_end) ?
      SelfIntersectionKind::kRepeatedVertex :
      SelfIntersectionKind::kNonAdjacentTouch;
  }

  if (first_first * first_second < 0 &&
      second_first * second_second < 0)
  {
    return SelfIntersectionKind::kProperCrossing;
  }

  const bool touches =
    (first_first == 0 && onClosedSegment(
      first_start, first_end, second_start)) ||
    (first_second == 0 && onClosedSegment(
      first_start, first_end, second_end)) ||
    (second_first == 0 && onClosedSegment(
      second_start, second_end, first_start)) ||
    (second_second == 0 && onClosedSegment(
      second_start, second_end, first_end));
  if (!touches)
    return SelfIntersectionKind::kNone;
  return segmentsShareEndpoint(
    first_start, first_end, second_start, second_end) ?
    SelfIntersectionKind::kRepeatedVertex :
    SelfIntersectionKind::kNonAdjacentTouch;
}

inline const char * selfIntersectionKindName(SelfIntersectionKind kind)
{
  switch (kind)
  {
    case SelfIntersectionKind::kProperCrossing:
      return "proper crossing";
    case SelfIntersectionKind::kNonAdjacentTouch:
      return "non-adjacent touch";
    case SelfIntersectionKind::kRepeatedVertex:
      return "repeated vertex";
    case SelfIntersectionKind::kCollinearOverlap:
      return "collinear overlap";
    case SelfIntersectionKind::kNone:
      return "no intersection";
  }
  return "unknown intersection";
}

inline SelfIntersectionResult selfIntersectionInconclusive(
  std::string diagnostic)
{
  SelfIntersectionResult result;
  result.status = SelfIntersectionStatus::kInconclusive;
  result.diagnostic = std::move(diagnostic);
  return result;
}

struct GridPoint
{
  long double x = 0.0L;
  long double y = 0.0L;
};

struct ParameterInterval
{
  long double lower = 0.0L;
  long double upper = 0.0L;
};

inline ValidationResult inconclusive(std::string diagnostic)
{
  ValidationResult result;
  result.status = ValidationStatus::kInconclusive;
  result.diagnostic = std::move(diagnostic);
  return result;
}

inline std::string formatGridPoint(const GridPoint & point)
{
  std::ostringstream stream;
  stream << std::setprecision(std::numeric_limits<long double>::max_digits10)
         << '(' << point.x << ", " << point.y << ')';
  return stream.str();
}

inline long double parameterTolerance(long double first, long double second)
{
  const long double scale = std::max(
    {1.0L, std::abs(first), std::abs(second)});
  return 64.0L * std::numeric_limits<long double>::epsilon() * scale;
}

// Recover an integer grid boundary only when the supplied double-precision
// world coordinate is within two double ULPs of that boundary's independent
// forward transform.  The actual world-to-grid division remains long double.
inline double boundarySnapTolerance(double boundary)
{
  if (!std::isfinite(boundary))
    return std::numeric_limits<double>::infinity();
  const double upward_ulp = std::abs(
    std::nextafter(boundary, std::numeric_limits<double>::infinity()) -
    boundary);
  const double downward_ulp = std::abs(
    boundary - std::nextafter(
      boundary, -std::numeric_limits<double>::infinity()));
  return 2.0 * std::max(upward_ulp, downward_ulp);
}

inline long double canonicalizeGridBoundary(
  double world, double origin, double resolution, std::uint32_t extent,
  long double converted)
{
  if (!std::isfinite(converted))
    return converted;

  const long double nearest = std::round(converted);
  if (!std::isfinite(nearest) || nearest < 0.0L ||
      nearest > static_cast<long double>(extent))
  {
    return converted;
  }

  const double reconstructed =
    origin + static_cast<double>(nearest) * resolution;
  if (!std::isfinite(reconstructed))
    return converted;

  const double tolerance = boundarySnapTolerance(reconstructed);
  if (world == reconstructed ||
      (std::isfinite(tolerance) &&
       std::abs(world - reconstructed) <= tolerance))
  {
    return nearest;
  }
  return converted;
}

inline bool hasSeparatedWorldBoundarySnapWindows(
  double origin, double resolution, std::uint32_t extent)
{
  double previous = origin;
  double previous_tolerance = boundarySnapTolerance(previous);
  if (!std::isfinite(previous_tolerance))
    return false;
  for (std::uint64_t index = 1; index <= extent; ++index)
  {
    const double current = origin + static_cast<double>(index) * resolution;
    const double current_tolerance = boundarySnapTolerance(current);
    const long double separation =
      static_cast<long double>(current) -
      static_cast<long double>(previous);
    const long double combined_tolerance =
      static_cast<long double>(previous_tolerance) +
      static_cast<long double>(current_tolerance);
    if (!std::isfinite(current) || !std::isfinite(current_tolerance) ||
        !(separation > combined_tolerance))
    {
      return false;
    }
    previous = current;
    previous_tolerance = current_tolerance;
  }
  return true;
}

inline ValidationResult validateMap(
  const nav_msgs::msg::OccupancyGrid & map, const OracleOptions & options)
{
  if (map.header.frame_id.empty())
    return inconclusive("Collision oracle requires a non-empty map frame_id");

  if (map.info.width == 0 || map.info.height == 0)
    return inconclusive("Collision oracle requires a non-empty map");

  const double resolution = static_cast<double>(map.info.resolution);
  if (!std::isfinite(resolution) || resolution <= 0.0)
    return inconclusive("Collision oracle requires a finite positive resolution");

  const auto & origin = map.info.origin;
  if (!std::isfinite(origin.position.x) ||
      !std::isfinite(origin.position.y) ||
      !std::isfinite(origin.position.z))
  {
    return inconclusive("Collision oracle requires a finite map origin");
  }
  if (origin.position.z != 0.0)
    return inconclusive("Collision oracle supports only a zero origin z coordinate");

  const auto & orientation = origin.orientation;
  if (!std::isfinite(orientation.x) || !std::isfinite(orientation.y) ||
      !std::isfinite(orientation.z) || !std::isfinite(orientation.w))
  {
    return inconclusive("Collision oracle requires a finite origin orientation");
  }
  const long double norm_squared =
    static_cast<long double>(orientation.x) * orientation.x +
    static_cast<long double>(orientation.y) * orientation.y +
    static_cast<long double>(orientation.z) * orientation.z +
    static_cast<long double>(orientation.w) * orientation.w;
  if (orientation.x != 0.0 || orientation.y != 0.0 ||
      orientation.z != 0.0 ||
      std::abs(norm_squared - 1.0L) > 1.0e-12L)
  {
    return inconclusive(
      "Collision oracle supports only an identity map rotation");
  }

  if (options.occupied_threshold < 1 || options.occupied_threshold > 100)
  {
    return inconclusive("occupied_threshold must be in [1, 100]");
  }

  const std::size_t width = static_cast<std::size_t>(map.info.width);
  const std::size_t height = static_cast<std::size_t>(map.info.height);
  if (width > std::numeric_limits<std::size_t>::max() / height)
    return inconclusive("OccupancyGrid dimensions overflow size_t");
  const std::size_t expected_cells = width * height;
  if (map.data.size() != expected_cells)
  {
    return inconclusive(
      "OccupancyGrid data length does not match width * height");
  }

  for (std::size_t index = 0; index < map.data.size(); ++index)
  {
    const int value = static_cast<int>(map.data[index]);
    if (value < -1 || value > 100)
    {
      const std::size_t x = index % width;
      const std::size_t y = index / width;
      return inconclusive(
        "Unsupported occupancy value " + std::to_string(value) +
        " at cell (" + std::to_string(x) + ", " + std::to_string(y) + ')');
    }
  }

  // If two adjacent grid boundaries collapse to the same double, a returned
  // ROS pose cannot identify which boundary it represents.  Report this as
  // inconclusive instead of guessing with a broad geometric tolerance.
  if (!hasSeparatedWorldBoundarySnapWindows(
      origin.position.x, resolution, map.info.width) ||
      !hasSeparatedWorldBoundarySnapWindows(
      origin.position.y, resolution, map.info.height))
  {
    return inconclusive(
      "Adjacent grid boundaries are not distinguishable outside their "
      "two-ULP snap windows");
  }

  ValidationResult result;
  result.status = ValidationStatus::kCollisionFree;
  result.diagnostic = "Map is suitable for independent collision validation";
  return result;
}

inline bool worldToGrid(
  const nav_msgs::msg::OccupancyGrid & map,
  const geometry_msgs::msg::Point & world, GridPoint & grid,
  std::string & error)
{
  if (!std::isfinite(world.x) || !std::isfinite(world.y) ||
      !std::isfinite(world.z))
  {
    error = "Path contains a non-finite pose position";
    return false;
  }
  if (world.z != 0.0)
  {
    error = "Collision oracle supports only path poses with z = 0";
    return false;
  }

  const double resolution = static_cast<double>(map.info.resolution);
  const long double long_resolution =
    static_cast<long double>(resolution);
  long double x =
    (static_cast<long double>(world.x) -
     static_cast<long double>(map.info.origin.position.x)) /
    long_resolution;
  long double y =
    (static_cast<long double>(world.y) -
     static_cast<long double>(map.info.origin.position.y)) /
    long_resolution;
  x = canonicalizeGridBoundary(
    world.x, map.info.origin.position.x, resolution, map.info.width, x);
  y = canonicalizeGridBoundary(
    world.y, map.info.origin.position.y, resolution, map.info.height, y);
  if (!std::isfinite(x) || !std::isfinite(y))
  {
    error = "World-to-grid conversion produced a non-finite coordinate";
    return false;
  }
  grid = GridPoint{x, y};
  error.clear();
  return true;
}

inline bool isFreeCell(
  const nav_msgs::msg::OccupancyGrid & map, std::uint32_t x,
  std::uint32_t y, const OracleOptions & options)
{
  // The planner deliberately forces the complete outer ring occupied,
  // independently of the input data values.
  if (x == 0 || y == 0 || x + 1 == map.info.width ||
      y + 1 == map.info.height)
  {
    return false;
  }

  const std::size_t index =
    static_cast<std::size_t>(y) * map.info.width + x;
  const int value = static_cast<int>(map.data[index]);
  if (value == -1)
    return options.allow_unknown;
  return value < options.occupied_threshold;
}

inline bool intersectClosedAxis(
  long double start, long double direction, long double minimum,
  long double maximum, long double & lower, long double & upper)
{
  if (direction == 0.0L)
    return start >= minimum && start <= maximum;

  long double first = (minimum - start) / direction;
  long double second = (maximum - start) / direction;
  if (first > second)
    std::swap(first, second);
  lower = std::max(lower, first);
  upper = std::min(upper, second);
  return true;
}

inline bool segmentCellClosureInterval(
  const GridPoint & start, const GridPoint & end, std::uint32_t cell_x,
  std::uint32_t cell_y, ParameterInterval & interval)
{
  long double lower = 0.0L;
  long double upper = 1.0L;
  const long double dx = end.x - start.x;
  const long double dy = end.y - start.y;
  if (!intersectClosedAxis(
      start.x, dx, static_cast<long double>(cell_x),
      static_cast<long double>(cell_x) + 1.0L, lower, upper) ||
      !intersectClosedAxis(
      start.y, dy, static_cast<long double>(cell_y),
      static_cast<long double>(cell_y) + 1.0L, lower, upper))
  {
    return false;
  }

  const long double tolerance = parameterTolerance(lower, upper);
  if (lower > upper + tolerance)
    return false;
  if (lower > upper)
  {
    const long double contact = std::clamp(
      (lower + upper) / 2.0L, 0.0L, 1.0L);
    lower = contact;
    upper = contact;
  }
  interval = ParameterInterval{
    std::clamp(lower, 0.0L, 1.0L),
    std::clamp(upper, 0.0L, 1.0L)};
  return true;
}

inline ValidationResult collisionResult(
  const GridPoint & start, const GridPoint & end, long double parameter,
  std::size_t segment_index, const std::string & reason)
{
  const long double clamped = std::clamp(parameter, 0.0L, 1.0L);
  const GridPoint collision{
    start.x + clamped * (end.x - start.x),
    start.y + clamped * (end.y - start.y)};

  ValidationResult result;
  result.status = ValidationStatus::kCollision;
  result.segment_index = segment_index;
  result.grid_x = collision.x;
  result.grid_y = collision.y;
  result.diagnostic = reason + " near grid coordinate " +
    formatGridPoint(collision);
  return result;
}

inline ValidationResult validateGridSegment(
  const nav_msgs::msg::OccupancyGrid & map, const GridPoint & start,
  const GridPoint & end, const OracleOptions & options,
  std::size_t segment_index)
{
  std::vector<ParameterInterval> free_intervals;
  free_intervals.reserve(map.data.size());
  for (std::uint32_t y = 0; y < map.info.height; ++y)
  {
    for (std::uint32_t x = 0; x < map.info.width; ++x)
    {
      if (!isFreeCell(map, x, y, options))
        continue;
      ParameterInterval interval;
      if (segmentCellClosureInterval(start, end, x, y, interval))
        free_intervals.push_back(interval);
    }
  }

  std::sort(
    free_intervals.begin(), free_intervals.end(),
    [](const ParameterInterval & lhs, const ParameterInterval & rhs) {
      if (lhs.lower != rhs.lower)
        return lhs.lower < rhs.lower;
      return lhs.upper > rhs.upper;
    });

  long double covered_until = 0.0L;
  bool coverage_started = false;
  for (const ParameterInterval & interval : free_intervals)
  {
    if (!coverage_started)
    {
      const long double tolerance = parameterTolerance(0.0L, interval.lower);
      if (interval.lower > tolerance)
      {
        return collisionResult(
          start, end, 0.0L, segment_index,
          "Segment starts outside the closure of every free cell");
      }
      coverage_started = true;
      covered_until = interval.upper;
    }
    else
    {
      const long double tolerance =
        parameterTolerance(covered_until, interval.lower);
      if (interval.lower > covered_until + tolerance)
      {
        return collisionResult(
          start, end, (covered_until + interval.lower) / 2.0L,
          segment_index,
          "Segment contains a gap not covered by any free-cell closure");
      }
      covered_until = std::max(covered_until, interval.upper);
    }

    if (covered_until >=
      1.0L - parameterTolerance(covered_until, 1.0L))
    {
      ValidationResult result;
      result.status = ValidationStatus::kCollisionFree;
      result.segment_index = segment_index;
      result.diagnostic =
        "Segment is completely covered by free-cell closures";
      return result;
    }
  }

  if (!coverage_started)
  {
    return collisionResult(
      start, end, 0.0L, segment_index,
      "Segment does not intersect the closure of any free cell");
  }
  return collisionResult(
    start, end, (covered_until + 1.0L) / 2.0L, segment_index,
    "Segment ends outside the closure of every free cell");
}

}  // namespace detail

inline ValidationResult validateSegment(
  const nav_msgs::msg::OccupancyGrid & map,
  const geometry_msgs::msg::Point & start,
  const geometry_msgs::msg::Point & end,
  const OracleOptions & options = OracleOptions{},
  std::size_t segment_index = 0)
{
  const ValidationResult map_result = detail::validateMap(map, options);
  if (!map_result.collisionFree())
    return map_result;

  detail::GridPoint grid_start;
  detail::GridPoint grid_end;
  std::string error;
  if (!detail::worldToGrid(map, start, grid_start, error) ||
      !detail::worldToGrid(map, end, grid_end, error))
  {
    ValidationResult result = detail::inconclusive(error);
    result.segment_index = segment_index;
    return result;
  }
  return detail::validateGridSegment(
    map, grid_start, grid_end, options, segment_index);
}

inline ValidationResult validatePath(
  const nav_msgs::msg::OccupancyGrid & map,
  const nav_msgs::msg::Path & path,
  const OracleOptions & options = OracleOptions{})
{
  const ValidationResult map_result = detail::validateMap(map, options);
  if (!map_result.collisionFree())
    return map_result;
  if (path.header.frame_id != map.header.frame_id)
  {
    return detail::inconclusive(
      "Path frame_id '" + path.header.frame_id +
      "' does not match map frame_id '" + map.header.frame_id + "'");
  }
  if (path.poses.empty())
    return detail::inconclusive("Cannot validate an empty path");

  for (std::size_t index = 0; index < path.poses.size(); ++index)
  {
    if (path.poses[index].header.frame_id != map.header.frame_id)
    {
      return detail::inconclusive(
        "Path pose " + std::to_string(index) + " frame_id '" +
        path.poses[index].header.frame_id +
        "' does not match map frame_id '" + map.header.frame_id + "'");
    }
  }

  if (path.poses.size() == 1)
  {
    detail::GridPoint point;
    std::string error;
    if (!detail::worldToGrid(map, path.poses.front().pose.position, point, error))
      return detail::inconclusive(error);
    return detail::validateGridSegment(map, point, point, options, 0);
  }

  for (std::size_t index = 0; index + 1 < path.poses.size(); ++index)
  {
    detail::GridPoint start;
    detail::GridPoint end;
    std::string error;
    if (!detail::worldToGrid(
        map, path.poses[index].pose.position, start, error) ||
        !detail::worldToGrid(
        map, path.poses[index + 1].pose.position, end, error))
    {
      ValidationResult result = detail::inconclusive(error);
      result.segment_index = index;
      return result;
    }
    ValidationResult result =
      detail::validateGridSegment(map, start, end, options, index);
    if (!result.collisionFree())
      return result;
  }

  ValidationResult result;
  result.status = ValidationStatus::kCollisionFree;
  result.segment_index = path.poses.size() - 2;
  result.diagnostic = "Every path segment passed the independent collision oracle";
  return result;
}

// Exhaustively inspect the final ROS polyline rather than the Core waypoint
// chain.  This catches regressions introduced by path reconstruction,
// interpolation, or ROS response assembly as well as planner pruning bugs.
inline SelfIntersectionResult validateNoSelfIntersection(
  const nav_msgs::msg::Path & path,
  std::size_t max_segments = 4096)
{
  if (path.header.frame_id.empty())
  {
    return detail::selfIntersectionInconclusive(
      "Self-intersection oracle requires a non-empty path frame_id");
  }
  if (path.poses.empty())
  {
    return detail::selfIntersectionInconclusive(
      "Cannot validate self-intersection for an empty path");
  }
  if (path.poses.size() - 1 > max_segments)
  {
    return detail::selfIntersectionInconclusive(
      "Self-intersection oracle segment count exceeds max_segments=" +
      std::to_string(max_segments));
  }

  std::vector<detail::ExactPlanarPoint> points;
  points.reserve(path.poses.size());
  for (std::size_t index = 0; index < path.poses.size(); ++index)
  {
    const auto & pose = path.poses[index];
    if (pose.header.frame_id != path.header.frame_id)
    {
      return detail::selfIntersectionInconclusive(
        "Path pose " + std::to_string(index) + " frame_id '" +
        pose.header.frame_id + "' does not match path frame_id '" +
        path.header.frame_id + "'");
    }
    const auto & position = pose.pose.position;
    if (!std::isfinite(position.x) || !std::isfinite(position.y) ||
        !std::isfinite(position.z))
    {
      return detail::selfIntersectionInconclusive(
        "Path pose " + std::to_string(index) +
        " contains a non-finite position");
    }
    if (position.z != 0.0)
    {
      return detail::selfIntersectionInconclusive(
        "Self-intersection oracle supports only path poses with z = 0");
    }
    points.emplace_back(
      detail::makeExactPlanarPoint(position.x, position.y));
  }

  for (std::size_t first = 0; first + 1 < points.size(); ++first)
  {
    // second = first + 1 is the consecutive segment.  Its shared construction
    // waypoint is legal and mirrors the production pruning contract.
    for (std::size_t second = first + 2;
         second + 1 < points.size(); ++second)
    {
      const SelfIntersectionKind kind =
        detail::classifyClosedSegmentIntersection(
          points[first], points[first + 1],
          points[second], points[second + 1]);
      if (kind == SelfIntersectionKind::kNone)
        continue;

      SelfIntersectionResult result;
      result.status = SelfIntersectionStatus::kSelfIntersection;
      result.kind = kind;
      result.first_segment_index = first;
      result.second_segment_index = second;
      result.diagnostic = std::string("Path has a ") +
        detail::selfIntersectionKindName(kind) + " between segments " +
        std::to_string(first) + " and " + std::to_string(second);
      return result;
    }
  }

  SelfIntersectionResult result;
  result.status = SelfIntersectionStatus::kIntersectionFree;
  result.kind = SelfIntersectionKind::kNone;
  result.diagnostic =
    "No pair of non-consecutive path segments intersects";
  return result;
}

}  // namespace raystar::test_oracle

#endif  // RAYSTAR_TEST_PATH_COLLISION_ORACLE_H_
