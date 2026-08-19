#ifndef RAYSTAR_POLYMAP_DETAIL_H_
#define RAYSTAR_POLYMAP_DETAIL_H_

// Package-private helpers shared by the Polymap translation units
// (polymap.cpp, polymap_visibility.cpp, polymap_shortening.cpp).
// Formerly an anonymous namespace inside the single polymap.cpp.

#include <raystar/polymap.h>
#include "visibility_validation.h"
#include <unordered_set>
#include <stack>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <functional>
#include <list>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <CGAL/exceptions.h>
#include <CGAL/number_utils.h>


namespace raystar {

namespace polymap_impl {


// Polymap is also a public direct-Core entry point, so it cannot rely solely
// on RaystarCore's request validator.  Validate dimensions and the backing
// vector before copying the map or sizing the two cell-sized topology
// registries.  The returned cell_count is always safe for the int-based
// indexing used throughout the legacy geometry code.
inline bool prepareGridStorage(
  const GridMap& grid_map, int& xsize, int& ysize, size_t& cell_count, std::string& error) {
  const size_t width = static_cast<size_t>(grid_map.width);
  const size_t height = static_cast<size_t>(grid_map.height);
  if (width == 0 || height == 0) {
    error = "Invalid map: width and height must be positive";
    return false;
  }
  if (width > std::numeric_limits<size_t>::max() / height) {
    error = "Invalid map: width * height overflows size_t";
    return false;
  }

  const size_t max_int = static_cast<size_t>(std::numeric_limits<int>::max());
  if (width > max_int || height > max_int) {
    error = "Invalid map: width and height must fit signed int indexing";
    return false;
  }
  cell_count = width * height;
  if (cell_count > max_int) {
    error = "Invalid map: cell count must fit signed int indexing";
    return false;
  }
  if (grid_map.data.size() != cell_count) {
    error = "Invalid map: data size is " + std::to_string(grid_map.data.size()) + ", expected " +
            std::to_string(cell_count);
    return false;
  }

  xsize = static_cast<int>(width);
  ysize = static_cast<int>(height);
  error.clear();
  return true;
}

inline bool validateCDTWithCGAL(const constrained_delaunay_triangulation::CDT& cdt) {
  return cdt.is_valid(false, 0);
}

struct UndirectedGridEdgeKey {
  int low;
  int high;

  bool operator==(const UndirectedGridEdgeKey& other) const noexcept {
    return low == other.low && high == other.high;
  }
};

inline UndirectedGridEdgeKey makeUndirectedGridEdgeKey(int first, int second) noexcept {
  if (first <= second)
    return {first, second};
  return {second, first};
}

struct UndirectedGridEdgeKeyHash {
  size_t operator()(const UndirectedGridEdgeKey& edge) const noexcept {
    const uint64_t packed = (static_cast<uint64_t>(static_cast<uint32_t>(edge.low)) << 32) |
                            static_cast<uint32_t>(edge.high);
    return std::hash<uint64_t>{}(packed);
  }
};

struct DirectedGridEdge {
  int from;
  int to;
};

inline long double cross2(long double ax, long double ay, long double bx, long double by) {
  return ax * by - ay * bx;
}

inline long double cross2(const Point2d& origin, const Point2d& first, const Point2d& second) {
  return cross2(static_cast<long double>(first.first) - origin.first,
                static_cast<long double>(first.second) - origin.second,
                static_cast<long double>(second.first) - origin.first,
                static_cast<long double>(second.second) - origin.second);
}

inline bool samePoint(const Point2d& first, const Point2d& second) {
  return first.first == second.first && first.second == second.second;
}

inline double polylineLength(const std::vector<Point2d>& path) {
  double length = 0.0;
  for (size_t index = 1; index < path.size(); ++index) {
    length += std::hypot(path[index].first - path[index - 1].first,
                         path[index].second - path[index - 1].second);
  }
  return length;
}

inline OperationStatus pointInPolygon(const std::vector<std::pair<int, int>>& polygon,
                               long double x,
                               long double y,
                               const StopToken& stop_token,
                               bool& inside) {
  inside = false;
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (polygon.size() < 3)
    return OperationStatus::success;
  for (size_t current = 0, previous = polygon.size() - 1; current < polygon.size();
       previous = current++) {
    if (stop_token.poll()) {
      inside = false;
      return OperationStatus::stopped;
    }
    const long double current_x = polygon[current].first;
    const long double current_y = polygon[current].second;
    const long double previous_x = polygon[previous].first;
    const long double previous_y = polygon[previous].second;
    const bool straddles = (current_y > y) != (previous_y > y);
    if (!straddles)
      continue;
    const long double intersection_x =
      current_x + (previous_x - current_x) * (y - current_y) / (previous_y - current_y);
    if (x < intersection_x)
      inside = !inside;
  }
  return OperationStatus::success;
}

// Raw reference contours coexist with the constrained triangulation, its
// directed-facet lookup, the edge accumulator, and the retained triangle
// mesh.  Their node-based CGAL/std containers have implementation-dependent
// allocator overhead, so max_map_bytes is deliberately charged a coarse 4 KiB
// admission estimate per unsimplified contour vertex.  This is a fail-closed
// complexity budget, not a claim about exact heap consumption.
constexpr size_t kEstimatedReferenceTopologyBytesPerRawContourVertex = 4096u;

inline bool pointInTriangleClosed(const std::array<Point2d, 3>& vertices, const Point2d& point) {
  const long double first = cross2(vertices[0], vertices[1], point);
  const long double second = cross2(vertices[1], vertices[2], point);
  const long double third = cross2(vertices[2], vertices[0], point);
  long double scale = 1.0L;
  for (const auto& vertex : vertices) {
    scale = std::max(scale,
                     std::abs(static_cast<long double>(vertex.first)) +
                       std::abs(static_cast<long double>(vertex.second)));
  }
  const long double tolerance = 1.0e-14L * scale * scale;
  const bool has_negative = first < -tolerance || second < -tolerance || third < -tolerance;
  const bool has_positive = first > tolerance || second > tolerance || third > tolerance;
  return !(has_negative && has_positive);
}

inline void appendSegmentEdgeEvents(const Point2d& start,
                             const Point2d& goal,
                             const Point2d& edge_start,
                             const Point2d& edge_goal,
                             std::vector<long double>& events) {
  const long double dx = static_cast<long double>(goal.first) - start.first;
  const long double dy = static_cast<long double>(goal.second) - start.second;
  const long double ex = static_cast<long double>(edge_goal.first) - edge_start.first;
  const long double ey = static_cast<long double>(edge_goal.second) - edge_start.second;
  const long double offset_x = static_cast<long double>(edge_start.first) - start.first;
  const long double offset_y = static_cast<long double>(edge_start.second) - start.second;
  const long double denominator = cross2(dx, dy, ex, ey);
  const long double scale =
    std::max<long double>(1.0L, std::max(std::abs(dx) + std::abs(dy), std::abs(ex) + std::abs(ey)));
  const long double tolerance = 1.0e-18L * scale * scale;

  const auto append_if_in_segment = [&events](long double parameter) {
    constexpr long double parameter_tolerance = 1.0e-14L;
    if (parameter < -parameter_tolerance || parameter > 1.0L + parameter_tolerance)
      return;
    events.push_back(std::max(0.0L, std::min(1.0L, parameter)));
  };

  if (std::abs(denominator) > tolerance) {
    const long double parameter = cross2(offset_x, offset_y, ex, ey) / denominator;
    const long double edge_parameter = cross2(offset_x, offset_y, dx, dy) / denominator;
    constexpr long double parameter_tolerance = 1.0e-14L;
    if (edge_parameter >= -parameter_tolerance && edge_parameter <= 1.0L + parameter_tolerance) {
      append_if_in_segment(parameter);
    }
    return;
  }

  if (std::abs(cross2(offset_x, offset_y, dx, dy)) > tolerance)
    return;
  const long double squared_length = dx * dx + dy * dy;
  if (squared_length == 0.0L)
    return;
  append_if_in_segment((offset_x * dx + offset_y * dy) / squared_length);
  const long double second_offset_x = static_cast<long double>(edge_goal.first) - start.first;
  const long double second_offset_y = static_cast<long double>(edge_goal.second) - start.second;
  append_if_in_segment((second_offset_x * dx + second_offset_y * dy) / squared_length);
}

inline int exactOrientation(const Point2d& first, const Point2d& second, const Point2d& third) {
  const auto orientation = CGAL::orientation(exact_geometry::Point(first.first, first.second),
                                             exact_geometry::Point(second.first, second.second),
                                             exact_geometry::Point(third.first, third.second));
  if (orientation == CGAL::LEFT_TURN)
    return 1;
  if (orientation == CGAL::RIGHT_TURN)
    return -1;
  return 0;
}

inline std::vector<Point2d> runFunnel(const Point2d& start,
                               const Point2d& goal,
                               const std::vector<DirectedPortal>& corridor_portals) {
  struct FunnelPortal {
    Point2d left;
    Point2d right;
  };
  std::vector<FunnelPortal> portals;
  portals.reserve(corridor_portals.size() + 2);
  portals.push_back({start, start});
  for (const auto& portal : corridor_portals) portals.push_back({portal.left, portal.right});
  portals.push_back({goal, goal});

  std::vector<Point2d> output;
  output.reserve(portals.size());
  output.push_back(start);
  Point2d apex = start;
  Point2d left = start;
  Point2d right = start;
  size_t apex_index = 0;
  size_t left_index = 0;
  size_t right_index = 0;

  for (size_t index = 1; index < portals.size(); ++index) {
    const Point2d new_left = portals[index].left;
    const Point2d new_right = portals[index].right;

    // exactOrientation uses the conventional positive-counter-clockwise
    // sign. Tightening the left leg moves clockwise; tightening the right
    // leg moves counter-clockwise.
    if (exactOrientation(apex, left, new_left) <= 0) {
      if (samePoint(apex, left) || exactOrientation(apex, right, new_left) > 0) {
        left = new_left;
        left_index = index;
      } else {
        if (!samePoint(output.back(), right))
          output.push_back(right);
        apex = right;
        apex_index = right_index;
        left = apex;
        right = apex;
        left_index = apex_index;
        right_index = apex_index;
        index = apex_index;
        continue;
      }
    }

    if (exactOrientation(apex, right, new_right) >= 0) {
      if (samePoint(apex, right) || exactOrientation(apex, left, new_right) < 0) {
        right = new_right;
        right_index = index;
      } else {
        if (!samePoint(output.back(), left))
          output.push_back(left);
        apex = left;
        apex_index = left_index;
        left = apex;
        right = apex;
        left_index = apex_index;
        right_index = apex_index;
        index = apex_index;
      }
    }
  }

  if (!samePoint(output.back(), goal))
    output.push_back(goal);
  return output;
}

enum class ExactSegmentRelation { disjoint, endpoint_touch, t_junction, proper_crossing, overlap };

struct ExactSegmentRelationResult {
  ExactSegmentRelation relation = ExactSegmentRelation::disjoint;
  std::optional<exact_geometry::Point> contact;
};

inline bool oppositeOrientations(CGAL::Orientation first, CGAL::Orientation second) {
  return (first == CGAL::LEFT_TURN && second == CGAL::RIGHT_TURN) ||
         (first == CGAL::RIGHT_TURN && second == CGAL::LEFT_TURN);
}

inline ExactSegmentRelationResult classifyExactSegments(const exact_geometry::Point& first_from,
                                                 const exact_geometry::Point& first_to,
                                                 const exact_geometry::Point& second_from,
                                                 const exact_geometry::Point& second_to) {
  const auto first_segment = exact_geometry::Kernel::Segment_2(first_from, first_to);
  const auto second_segment = exact_geometry::Kernel::Segment_2(second_from, second_to);

  const CGAL::Orientation second_from_side = CGAL::orientation(first_from, first_to, second_from);
  const CGAL::Orientation second_to_side = CGAL::orientation(first_from, first_to, second_to);
  const CGAL::Orientation first_from_side = CGAL::orientation(second_from, second_to, first_from);
  const CGAL::Orientation first_to_side = CGAL::orientation(second_from, second_to, first_to);

  if (oppositeOrientations(second_from_side, second_to_side) &&
      oppositeOrientations(first_from_side, first_to_side)) {
    return {ExactSegmentRelation::proper_crossing, std::nullopt};
  }

  std::optional<exact_geometry::Point> common_endpoint;
  bool has_distinct_second_endpoint = false;
  const auto add_if_common = [&](const exact_geometry::Point& point) {
    if (!first_segment.has_on(point) || !second_segment.has_on(point))
      return;
    if (!common_endpoint)
      common_endpoint = point;
    else if (*common_endpoint != point)
      has_distinct_second_endpoint = true;
  };
  add_if_common(first_from);
  add_if_common(first_to);
  add_if_common(second_from);
  add_if_common(second_to);

  if (!common_endpoint)
    return {ExactSegmentRelation::disjoint, std::nullopt};
  if (has_distinct_second_endpoint)
    return {ExactSegmentRelation::overlap, std::nullopt};

  const auto& contact = *common_endpoint;
  const bool first_endpoint = contact == first_from || contact == first_to;
  const bool second_endpoint = contact == second_from || contact == second_to;
  if (first_endpoint && second_endpoint)
    return {ExactSegmentRelation::endpoint_touch, contact};
  return {ExactSegmentRelation::t_junction, contact};
}

inline const char* segmentRelationDescription(ExactSegmentRelation relation) {
  switch (relation) {
  case ExactSegmentRelation::endpoint_touch:
    return "share an endpoint";
  case ExactSegmentRelation::t_junction:
    return "form a T-junction";
  case ExactSegmentRelation::proper_crossing:
    return "properly cross";
  case ExactSegmentRelation::overlap:
    return "overlap";
  case ExactSegmentRelation::disjoint:
    break;
  }
  return "are disjoint";
}

std::optional<detail::VisibilityGeometryContext> makeVisibilityGeometryContext(
  const Polymap& polymap,
  const Point2d& source,
  const exact_geometry::Point& exact_source,
  const std::optional<std::pair<int, int>>& obstacle_vertex,
  std::string& error);

inline std::optional<detail::VisibilityGeometryContext> makeVisibilityGeometryContext(
  const Polymap& polymap, int source_x, int source_y, std::string& error) {
  const auto source_topology = polymap.locateVertex(source_x, source_y);
  return makeVisibilityGeometryContext(
    polymap,
    Point2d{static_cast<double>(source_x), static_cast<double>(source_y)},
    exact_geometry::Point(source_x, source_y),
    Polymap::isNoTopology(source_topology) ? std::optional<std::pair<int, int>>()
                                           : std::optional<std::pair<int, int>>(source_topology),
    error);
}

// Build the source-dependent validation context for both integer obstacle
// vertices and continuous free-space roots.  The optional topology is an
// explicit signal: a continuous point that happens to round to an obstacle
// vertex must still use the closed-cycle semantics (and is rejected by the
// endpoint validator before planning).
inline std::optional<detail::VisibilityGeometryContext> makeVisibilityGeometryContext(
  const Polymap& polymap,
  const Point2d& source,
  const exact_geometry::Point& exact_source,
  const std::optional<std::pair<int, int>>& obstacle_vertex,
  std::string& error) {
  error.clear();
  detail::VisibilityGeometryContext context{
    exact_source, detail::VisibilityBoundaryMode::closed_cycle, std::nullopt, std::nullopt};

  if (!std::isfinite(source.first) || !std::isfinite(source.second)) {
    error = "Visibility source is not finite";
    return std::nullopt;
  }

  if (!obstacle_vertex)
    return context;

  if (!polymap.isValidTopology(*obstacle_vertex)) {
    error = "Visibility source has an invalid obstacle topology";
    return std::nullopt;
  }

  const auto previous = polymap.getPrevObs(*obstacle_vertex);
  const auto next = polymap.getNextObs(*obstacle_vertex);
  if (!previous || !next) {
    error = "Visibility source obstacle has no valid sector anchors";
    return std::nullopt;
  }

  context.mode = detail::VisibilityBoundaryMode::open_sector;
  context.start_anchor.emplace(previous->first, previous->second);
  context.end_anchor.emplace(next->first, next->second);
  return context;
}

enum class ExactContourPointLocation { outside, boundary, inside };

inline ExactContourPointLocation locateExactPointInContour(
  const std::vector<std::pair<int, int>>& vertices,
  const exact_geometry::Point& point,
  const StopToken& stop_token) {
  if (vertices.size() < 3)
    return ExactContourPointLocation::outside;

  int winding_number = 0;
  for (size_t index = 0; index < vertices.size(); ++index) {
    if (stop_token.poll())
      return ExactContourPointLocation::outside;
    const auto& from_coordinate = vertices[index];
    const auto& to_coordinate = vertices[(index + 1) % vertices.size()];
    const exact_geometry::Point from(from_coordinate.first, from_coordinate.second);
    const exact_geometry::Point to(to_coordinate.first, to_coordinate.second);
    const exact_geometry::Kernel::Segment_2 edge(from, to);
    if (edge.has_on(point))
      return ExactContourPointLocation::boundary;

    // Exact winding-number test.  The half-open y intervals avoid counting a
    // vertex twice, while orientation supplies an epsilon-free x comparison.
    const bool upward = from.y() <= point.y() && point.y() < to.y();
    const bool downward = to.y() <= point.y() && point.y() < from.y();
    const CGAL::Orientation orientation = CGAL::orientation(from, to, point);
    if (upward && orientation == CGAL::LEFT_TURN)
      ++winding_number;
    else if (downward && orientation == CGAL::RIGHT_TURN)
      --winding_number;
  }
  return winding_number == 0 ? ExactContourPointLocation::outside
                             : ExactContourPointLocation::inside;
}

// A contour set extracted from a map without an occupied outer border can
// contain only obstacle holes.  In that legacy representation there is no
// finite outer contour against which a strict free-space point can be
// classified; callers such as RaystarCore add the border before construction.
// Only perform the early endpoint check when a clockwise outer contour is
// actually present.
inline bool hasClockwiseOuterContour(const std::vector<Obs>& obstacles, const StopToken& stop_token) {
  for (const auto& obstacle : obstacles) {
    if (stop_token.poll())
      return false;
    exact_geometry::FT twice_area(0);
    for (size_t index = 0; index < obstacle.ordered_vertices_.size(); ++index) {
      if (stop_token.poll())
        return false;
      const auto& from = obstacle.ordered_vertices_[index];
      const auto& to = obstacle.ordered_vertices_[(index + 1) % obstacle.ordered_vertices_.size()];
      twice_area +=
        exact_geometry::FT(from.first) * to.second - exact_geometry::FT(to.first) * from.second;
    }
    if (twice_area < exact_geometry::FT(0))
      return true;
  }
  return false;
}


}  // namespace polymap_impl

}  // namespace raystar

#endif  // RAYSTAR_POLYMAP_DETAIL_H_
