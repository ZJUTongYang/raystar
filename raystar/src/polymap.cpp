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

namespace {

// Polymap is also a public direct-Core entry point, so it cannot rely solely
// on RaystarCore's request validator.  Validate dimensions and the backing
// vector before copying the map or sizing the two cell-sized topology
// registries.  The returned cell_count is always safe for the int-based
// indexing used throughout the legacy geometry code.
bool prepareGridStorage(
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

bool validateCDTWithCGAL(const constrained_delaunay_triangulation::CDT& cdt) {
  return cdt.is_valid(false, 0);
}

struct UndirectedGridEdgeKey {
  int low;
  int high;

  bool operator==(const UndirectedGridEdgeKey& other) const noexcept {
    return low == other.low && high == other.high;
  }
};

UndirectedGridEdgeKey makeUndirectedGridEdgeKey(int first, int second) noexcept {
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

long double cross2(long double ax, long double ay, long double bx, long double by) {
  return ax * by - ay * bx;
}

long double cross2(const Point2d& origin, const Point2d& first, const Point2d& second) {
  return cross2(static_cast<long double>(first.first) - origin.first,
                static_cast<long double>(first.second) - origin.second,
                static_cast<long double>(second.first) - origin.first,
                static_cast<long double>(second.second) - origin.second);
}

bool samePoint(const Point2d& first, const Point2d& second) {
  return first.first == second.first && first.second == second.second;
}

double polylineLength(const std::vector<Point2d>& path) {
  double length = 0.0;
  for (size_t index = 1; index < path.size(); ++index) {
    length += std::hypot(path[index].first - path[index - 1].first,
                         path[index].second - path[index - 1].second);
  }
  return length;
}

OperationStatus pointInPolygon(const std::vector<std::pair<int, int>>& polygon,
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

bool pointInTriangleClosed(const std::array<Point2d, 3>& vertices, const Point2d& point) {
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

void appendSegmentEdgeEvents(const Point2d& start,
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

int exactOrientation(const Point2d& first, const Point2d& second, const Point2d& third) {
  const auto orientation = CGAL::orientation(exact_geometry::Point(first.first, first.second),
                                             exact_geometry::Point(second.first, second.second),
                                             exact_geometry::Point(third.first, third.second));
  if (orientation == CGAL::LEFT_TURN)
    return 1;
  if (orientation == CGAL::RIGHT_TURN)
    return -1;
  return 0;
}

std::vector<Point2d> runFunnel(const Point2d& start,
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

bool oppositeOrientations(CGAL::Orientation first, CGAL::Orientation second) {
  return (first == CGAL::LEFT_TURN && second == CGAL::RIGHT_TURN) ||
         (first == CGAL::RIGHT_TURN && second == CGAL::LEFT_TURN);
}

ExactSegmentRelationResult classifyExactSegments(const exact_geometry::Point& first_from,
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

const char* segmentRelationDescription(ExactSegmentRelation relation) {
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

std::optional<detail::VisibilityGeometryContext> makeVisibilityGeometryContext(
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
std::optional<detail::VisibilityGeometryContext> makeVisibilityGeometryContext(
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

ExactContourPointLocation locateExactPointInContour(
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
bool hasClockwiseOuterContour(const std::vector<Obs>& obstacles, const StopToken& stop_token) {
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

}  // namespace

bool validateReducedDirectedPortalWitness(const TriangleCorridor& corridor, std::string* error) {
  if (error)
    error->clear();
  const auto fail = [&](const std::string& message) {
    if (error)
      *error = message;
    return false;
  };

  if (corridor.triangle_occurrences.empty())
    return fail("A lifted-sleeve witness requires a containing triangle occurrence");
  if (corridor.portals.size() == std::numeric_limits<size_t>::max() ||
      corridor.triangle_occurrences.size() != corridor.portals.size() + 1) {
    return fail("Triangle-occurrence and directed-portal cardinalities are inconsistent");
  }

  for (size_t index = 0; index < corridor.portals.size(); ++index) {
    const auto& portal = corridor.portals[index];
    if (portal.from_triangle != corridor.triangle_occurrences[index] ||
        portal.to_triangle != corridor.triangle_occurrences[index + 1]) {
      return fail("Directed portal " + std::to_string(index) +
                  " does not bind its adjacent triangle occurrences");
    }
    if (portal.from_triangle == portal.to_triangle) {
      return fail("Directed portal " + std::to_string(index) + " is a triangle self-transition");
    }
    if (!std::isfinite(portal.left.first) || !std::isfinite(portal.left.second) ||
        !std::isfinite(portal.right.first) || !std::isfinite(portal.right.second)) {
      return fail("Directed portal " + std::to_string(index) + " contains a non-finite endpoint");
    }
    if (portal.left == portal.right) {
      return fail("Directed portal " + std::to_string(index) + " has zero geometric width");
    }
  }

  // traceFreeSpacePath() reduces the dual walk as a stack.  Retain later
  // repeated occurrences (which encode winding), but reject an adjacent
  // portal followed immediately by its inverse A->B->A.
  for (size_t index = 2; index < corridor.triangle_occurrences.size(); ++index) {
    if (corridor.triangle_occurrences[index] == corridor.triangle_occurrences[index - 2]) {
      return fail("Lifted-sleeve witness contains an unreduced immediate portal reversal at " +
                  std::to_string(index - 2));
    }
  }
  return true;
}

bool sameReducedDirectedPortalWitness(const TriangleCorridor& reference,
                                      const TriangleCorridor& candidate,
                                      std::string* error) {
  if (error)
    error->clear();
  const auto fail = [&](const std::string& message) {
    if (error)
      *error = message;
    return false;
  };

  std::string validation_error;
  if (!validateReducedDirectedPortalWitness(reference, &validation_error))
    return fail("Reference lifted-sleeve witness is invalid: " + validation_error);
  if (!validateReducedDirectedPortalWitness(candidate, &validation_error))
    return fail("Candidate lifted-sleeve witness is invalid: " + validation_error);

  if (reference.triangle_occurrences.size() != candidate.triangle_occurrences.size()) {
    return fail("Lifted-sleeve triangle-occurrence counts differ");
  }
  for (size_t index = 0; index < reference.triangle_occurrences.size(); ++index) {
    if (reference.triangle_occurrences[index] != candidate.triangle_occurrences[index]) {
      return fail("Lifted-sleeve triangle occurrence differs at ordinal " + std::to_string(index));
    }
  }

  // Cardinalities were validated above and the equal triangle counts imply
  // equal portal counts.  Compare each occurrence in order; never sort or
  // deduplicate repeated portal geometry.
  for (size_t index = 0; index < reference.portals.size(); ++index) {
    const auto& expected = reference.portals[index];
    const auto& actual = candidate.portals[index];
    if (expected.from_triangle != actual.from_triangle ||
        expected.to_triangle != actual.to_triangle) {
      return fail("Directed portal identity differs at occurrence " + std::to_string(index));
    }
    if (expected.left != actual.left || expected.right != actual.right) {
      return fail("Directed portal geometry differs at occurrence " + std::to_string(index));
    }
  }
  return true;
}

static std::optional<exact_geometry::Point> findIntersection(const exact_geometry::Point& s,
                                                             const exact_geometry::Point& g,
                                                             const exact_geometry::Point& p,
                                                             const exact_geometry::Point& limit) {
  return exact_geometry::intersectSegmentWithRay(s, g, p, limit);
}

static Point2d projectExactPoint(const exact_geometry::Point& point) {
  return {CGAL::to_double(point.x()), CGAL::to_double(point.y())};
}

// ============================================================================
// calculateVisibilityRegion
//
// Computes the visibility polygon from source point (round_x, round_y).
// Implements the Bungiu et al. 2014 algorithm using a CDT-based sweeper.
//
// Output:
//   visibility_region : polygon endpoints, with coordinates and obstacle
//                       support kept together as one value.
//
// BungiuEdge fields:
//   prev / next             : rich edge endpoints
//   limit_prev / limit_next : compatibility/debug angle projections
//   limit_*_anchor          : exact witnesses for the angular limits
//   is_bd                   : final boundary (skip expansion)
//   supporting_segment      : set only when the whole segment lies on one
//                             directed obstacle edge
//
// Main loop logic:
//   For each non-boundary edge (is_bd==false), find the adjacent CDT triangle.
//   The third vertex V is classified against the two directed limit rays with
//   EPECK orientation/dot-product predicates. Numeric angles are never used to
//   select a topology branch.
// ============================================================================
bool Polymap::calculateVisibilityRegion(int round_x,
                                        int round_y,
                                        VisibilityRegion& visibility_region) {
  const StopToken no_stop;
  return calculateVisibilityRegionImpl(round_x, round_y, visibility_region, no_stop);
}

OperationStatus Polymap::calculateVisibilityRegion(int round_x,
                                                   int round_y,
                                                   VisibilityRegion& visibility_region,
                                                   const StopToken& stop_token) {
  visibility_region.clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (calculateVisibilityRegionImpl(round_x, round_y, visibility_region, stop_token)) {
    if (stop_token.poll()) {
      visibility_region.clear();
      return OperationStatus::stopped;
    }
    return OperationStatus::success;
  }
  visibility_region.clear();
  return stop_token.poll() ? OperationStatus::stopped : OperationStatus::failure;
}

bool Polymap::calculateVisibilityRegionImpl(int round_x,
                                            int round_y,
                                            VisibilityRegion& visibility_region,
                                            const StopToken& stop_token) {
  const auto topology = locateVertex(round_x, round_y);
  const std::optional<std::pair<int, int>> obstacle_vertex =
    isValidTopology(topology) ? std::optional<std::pair<int, int>>(topology) : std::nullopt;
  return calculateVisibilityRegionImpl(
    Point2d{static_cast<double>(round_x), static_cast<double>(round_y)},
    exact_geometry::Point(round_x, round_y),
    obstacle_vertex,
    visibility_region,
    stop_token);
}

bool Polymap::calculateVisibilityRegionImpl(
  const Point2d& source,
  const exact_geometry::Point& exact_source,
  const std::optional<std::pair<int, int>>& obstacle_vertex_source,
  VisibilityRegion& visibility_region,
  const StopToken& stop_token) {
  visibility_region.clear();

  if (stop_token.poll())
    return false;

  if (!solution_exist_ || !cdt_ready_ || xsize_ <= 0 || ysize_ <= 0 ||
      !std::isfinite(source.first) || !std::isfinite(source.second))
    return false;

  const size_t expected_registry_size = static_cast<size_t>(xsize_) * static_cast<size_t>(ysize_);
  if (vertices_location_x_flat_.size() != expected_registry_size ||
      vertices_location_y_flat_.size() != expected_registry_size)
    return false;

  // The exact source is supplied by the caller so continuous root expansion
  // never goes through the lossy integer cache or a rounded coordinate.
  const double source_x = source.first;
  const double source_y = source.second;

  std::list<constrained_delaunay_triangulation::BungiuEdge> bd;
  bool open_visibility_region = false;
  std::pair<int, int> open_region_closing_vertex;
  std::pair<int, int> open_region_start_vertex;
  auto iter = bd.end();

  using BungiuEdge = constrained_delaunay_triangulation::BungiuEdge;
  const auto make_final_boundary = [](const BoundaryEndpoint& prev, const BoundaryEndpoint& next) {
    BungiuEdge edge;
    edge.prev = prev;
    edge.next = next;
    edge.is_bd = true;
    return edge;
  };
  const auto make_obstacle_boundary = [](const BoundaryEndpoint& prev,
                                         const BoundaryEndpoint& next,
                                         const DirectedObstacleEdge& support) {
    BungiuEdge edge;
    edge.prev = prev;
    edge.next = next;
    edge.is_bd = true;
    edge.supporting_segment = support;
    return edge;
  };
  const auto make_portal = [](const BoundaryEndpoint& prev,
                              const BoundaryEndpoint& next,
                              double limit_prev,
                              double limit_next,
                              const BoundaryEndpoint& limit_prev_anchor,
                              const BoundaryEndpoint& limit_next_anchor) {
    BungiuEdge edge;
    edge.prev = prev;
    edge.next = next;
    edge.limit_prev = limit_prev;
    edge.limit_next = limit_next;
    edge.limit_prev_anchor = limit_prev_anchor;
    edge.limit_next_anchor = limit_next_anchor;
    return edge;
  };
  const auto set_final_boundary =
    [](BungiuEdge& edge, const BoundaryEndpoint& prev, const BoundaryEndpoint& next) {
      edge.prev = prev;
      edge.next = next;
      edge.is_bd = true;
      edge.supporting_segment.reset();
    };
  const auto make_final_or_supported_boundary =
    [&](const BoundaryEndpoint& prev,
        const BoundaryEndpoint& next,
        const std::optional<DirectedObstacleEdge>& possible_support) {
      if (possible_support && pointLiesOnObstacleEdge(prev, *possible_support) &&
          pointLiesOnObstacleEdge(next, *possible_support)) {
        return make_obstacle_boundary(prev, next, *possible_support);
      }
      return make_final_boundary(prev, next);
    };
  const auto set_final_or_supported_boundary =
    [&](BungiuEdge& edge,
        const BoundaryEndpoint& prev,
        const BoundaryEndpoint& next,
        const std::optional<DirectedObstacleEdge>& possible_support) {
      edge = make_final_or_supported_boundary(prev, next, possible_support);
    };
  const auto validate_bungiu_edge = [&](const BungiuEdge& edge) {
    if (!validateBoundaryEndpointImpl(edge.prev, stop_token, nullptr) ||
        !validateBoundaryEndpointImpl(edge.next, stop_token, nullptr)) {
      return false;
    }
    if (edge.supporting_segment) {
      if (!edge.is_bd || !pointLiesOnObstacleEdge(edge.prev, *edge.supporting_segment) ||
          !pointLiesOnObstacleEdge(edge.next, *edge.supporting_segment)) {
        return false;
      }
    }
    if (!edge.is_bd) {
      if (edge.supporting_segment ||
          !validateBoundaryEndpointImpl(edge.limit_prev_anchor, stop_token, nullptr) ||
          !validateBoundaryEndpointImpl(edge.limit_next_anchor, stop_token, nullptr) ||
          !exact_geometry::isPortalSpanAtMostPi(
            exact_source, exactPoint(edge.limit_prev_anchor), exactPoint(edge.limit_next_anchor))) {
        return false;
      }
    }
    return true;
  };

  // --- Phase 1: Initialize boundary edges ---
  // Two cases: source IS an obstacle vertex, or source is in free space.

  if (obstacle_vertex_source) {
    // Case A: source is an obstacle vertex.
    // Walk around the source vertex through adjacent CDT triangles,
    // creating one BungiuEdge per triangle edge that radiates from the source.
    open_visibility_region = true;
    const auto curr = *obstacle_vertex_source;
    const auto source_coordinate =
      obs_[static_cast<size_t>(curr.first)].ordered_vertices_[static_cast<size_t>(curr.second)];
    // An obstacle source is always represented by its exact grid vertex.  Do
    // not silently combine a topology id with a different continuous point.
    if (exact_source != exact_geometry::Point(source_coordinate.first, source_coordinate.second))
      return false;
    auto prev_result = getPrevObs(curr);
    auto next_result = getNextObs(curr);
    if (!prev_result || !next_result)
      return false;
    const auto prev = *prev_result;
    const auto next = *next_result;
    open_region_closing_vertex = next;
    open_region_start_vertex = prev;

    auto the_vertex = prev;
    while (!(the_vertex.first == next.first && the_vertex.second == next.second)) {
      if (stop_token.poll())
        return false;
      int loc = locateAdjacentFacet(source_coordinate, the_vertex);
      if (loc < 0)
        return false;
      int i;
      for (i = 0; i < 3; ++i) {
        if (!(facets_[loc][i].first == source_coordinate.first &&
              facets_[loc][i].second == source_coordinate.second) &&
            !(facets_[loc][i].first == the_vertex.first &&
              facets_[loc][i].second == the_vertex.second))
          break;
      }
      if (i == 3)
        return false;
      bd.emplace_back(constrained_delaunay_triangulation::BungiuEdge());
      bd.back().prev = makeEndpoint(the_vertex);
      bd.back().next = makeEndpoint(facets_[loc][i]);
      bd.back().limit_prev =
        atan2(bd.back().prev.position.second - source_y, bd.back().prev.position.first - source_x);
      bd.back().limit_next =
        atan2(bd.back().next.position.second - source_y, bd.back().next.position.first - source_x);
      bd.back().limit_next = bd.back().limit_prev +
                             normalize_angle_positive(bd.back().limit_next - bd.back().limit_prev);
      bd.back().limit_prev_anchor = bd.back().prev;
      bd.back().limit_next_anchor = bd.back().next;
      bd.back().supporting_segment =
        obstacleEdgeBetween(bd.back().next.position, bd.back().prev.position);
      bd.back().is_bd = bd.back().supporting_segment.has_value();
      the_vertex = facets_[loc][i];
    }
  } else {
    // Case B: source is in free space (inside a CDT triangle).
    // Find the containing triangle and use its 3 edges as initial boundary.
    int loc = -1;
    for (size_t fi = 0; fi < facets_.size(); ++fi) {
      if (stop_token.poll())
        return false;
      if (isInTri(facets_[fi][0].first,
                  facets_[fi][0].second,
                  facets_[fi][1].first,
                  facets_[fi][1].second,
                  facets_[fi][2].first,
                  facets_[fi][2].second,
                  source_x,
                  source_y)) {
        loc = static_cast<int>(fi);
        break;
      }
    }
    if (loc < 0)
      return false;
    for (int i = 0; i < 3; ++i) {
      bd.emplace_back(constrained_delaunay_triangulation::BungiuEdge());
      bd.back().prev = makeEndpoint(facets_[loc][i]);
      bd.back().next = makeEndpoint(facets_[loc][(i + 1) % 3]);
      bd.back().limit_prev =
        atan2(bd.back().prev.position.second - source_y, bd.back().prev.position.first - source_x);
      bd.back().limit_next =
        atan2(bd.back().next.position.second - source_y, bd.back().next.position.first - source_x);
      bd.back().limit_next = bd.back().limit_prev +
                             normalize_angle_positive(bd.back().limit_next - bd.back().limit_prev);
      bd.back().limit_prev_anchor = bd.back().prev;
      bd.back().limit_next_anchor = bd.back().next;
      bd.back().supporting_segment =
        obstacleEdgeBetween(bd.back().next.position, bd.back().prev.position);
      bd.back().is_bd = bd.back().supporting_segment.has_value();
    }
  }

  // --- Phase 2: Expand boundary by sweeping through CDT triangles ---
  iter = bd.begin();
  double theta;
  exact_geometry::PortalRayPosition ray_position;
  int loc;
  while (1) {
    if (stop_token.poll())
      return false;
    if (iter == bd.end())
      break;

    if (!validate_bungiu_edge(*iter))
      return false;

    // Skip obstacle boundary edges - they are final, no triangle to expand into.
    if (iter->is_bd) {
      ++iter;
      continue;
    }

    // Find the CDT triangle adjacent to this edge (on the outer side).
    loc = locateAdjacentFacet(iter->next.position, iter->prev.position);
    if (loc < 0) {
      iter->is_bd = true;
      ++iter;
      continue;
    }
    int i;
    // Find the third vertex of the triangle (not prev_pos, not next_pos).
    for (i = 0; i < 3; ++i) {
      if (!(facets_[loc][i].first == iter->prev.position.first &&
            facets_[loc][i].second == iter->prev.position.second) &&
          !(facets_[loc][i].first == iter->next.position.first &&
            facets_[loc][i].second == iter->next.position.second))
        break;
    }
    if (i == 3)
      return false;
    const BoundaryEndpoint V = makeEndpoint(facets_[loc][i]);

    // Angle from source to the third vertex.
    theta = iter->limit_prev + normalize_angle(atan2(facets_[loc][i].second - source_y,
                                                     facets_[loc][i].first - source_x) -
                                               iter->limit_prev);
    ray_position = exact_geometry::classifyPortalRay(exact_source,
                                                     exactPoint(iter->limit_prev_anchor),
                                                     exactPoint(iter->limit_next_anchor),
                                                     exactPoint(V));
    if (ray_position == exact_geometry::PortalRayPosition::invalid)
      return false;

    if (ray_position == exact_geometry::PortalRayPosition::before_lower) {
      // CASE 1: theta < limit_prev
      // The triangle vertex V is at a smaller angle than prev_pos.
      // The obstacle edge from V to next_pos would block the view if they are consecutive.
      const auto blocker = obstacleEdgeBetween(iter->next.position, V.position);
      if (blocker) {
        // next_pos and V are connected by an obstacle edge.
        // The obstacle blocks the sweeper: compute far-endpoints where rays from
        // source through limit_prev_pos and limit_next_pos hit this obstacle edge.
        const auto endpoint_prev_position = findIntersection(
          exactPoint(V), exactPoint(iter->next), exact_source, exactPoint(iter->limit_prev_anchor));
        const auto endpoint_next_position = findIntersection(
          exactPoint(V), exactPoint(iter->next), exact_source, exactPoint(iter->limit_next_anchor));
        if (!endpoint_prev_position || !endpoint_next_position)
          return false;
        const auto endpoint_prev = makeIntersectionEndpoint(*endpoint_prev_position, *blocker);
        const auto endpoint_next = makeIntersectionEndpoint(*endpoint_next_position, *blocker);

        // Edge 1: limit_prev_pos -> endpoint_prev (boundary, ray hit point on prev side)
        bd.insert(
          iter, make_final_or_supported_boundary(iter->limit_prev_anchor, endpoint_prev, blocker));

        // Edge 2: endpoint_prev -> endpoint_next (along obstacle edge)
        bd.insert(iter, make_obstacle_boundary(endpoint_prev, endpoint_next, *blocker));

        // Update or erase the current edge
        if (exactPoint(iter->next) != exactPoint(endpoint_next)) {
          // endpoint_next is far from next_pos: keep iter as boundary edge
          set_final_or_supported_boundary(*iter, endpoint_next, iter->limit_next_anchor, blocker);
          iter++;
        } else {
          // endpoint_next coincides with next_pos: remove this edge
          bd.erase(iter++);
        }
      } else {
        // V lies outside the portal's lower angular limit. Advance the portal
        // endpoint through the adjacent triangle, but do not emit V or the
        // edge (V, prev_pos): that geometry is outside this angular window.
        iter->prev = V;
      }
    } else if (ray_position == exact_geometry::PortalRayPosition::equal_lower) {
      // CASE 2: theta == limit_prev
      // V is collinear with source and prev_pos (same angle).
      // This means prev_pos and V are at the same angle from source.
      const auto blocker = obstacleEdgeBetween(iter->next.position, V.position);
      // Insert boundary edge: limit_prev_pos -> V
      bd.insert(iter, make_final_or_supported_boundary(iter->limit_prev_anchor, V, blocker));
      if (blocker) {
        // Obstacle edge from next_pos to V: compute far-endpoint on next side.
        const auto endpoint_next_position = findIntersection(
          exactPoint(V), exactPoint(iter->next), exact_source, exactPoint(iter->limit_next_anchor));
        if (!endpoint_next_position)
          return false;
        const auto endpoint_next = makeIntersectionEndpoint(*endpoint_next_position, *blocker);
        bd.insert(iter, make_obstacle_boundary(V, endpoint_next, *blocker));

        if (exactPoint(iter->next) != exactPoint(endpoint_next)) {
          set_final_or_supported_boundary(*iter, endpoint_next, iter->limit_next_anchor, blocker);
          iter++;
        } else {
          bd.erase(iter++);
        }
      } else {
        // No obstacle edge: V extends the boundary, update prev.
        iter->prev = V;
        iter->limit_prev_anchor = V;
      }
    } else if (ray_position == exact_geometry::PortalRayPosition::inside) {
      // CASE 3: limit_prev < theta < limit_next
      // V is strictly inside the current angular range.
      // The current edge is split into two parts at V.
      bool to_minus = false;
      const auto prev_blocker = obstacleEdgeBetween(V.position, iter->prev.position);
      if (prev_blocker) {
        // Obstacle edge between V and prev_pos.
        // Compute far-endpoint where ray from source through limit_prev_pos hits this edge.
        const auto endpoint_prev_position = findIntersection(
          exactPoint(iter->prev), exactPoint(V), exact_source, exactPoint(iter->limit_prev_anchor));
        if (!endpoint_prev_position)
          return false;
        const auto endpoint_prev = makeIntersectionEndpoint(*endpoint_prev_position, *prev_blocker);

        if (exactPoint(endpoint_prev) != exactPoint(iter->limit_prev_anchor)) {
          // Far-endpoint is distinct from limit_prev_pos:
          // insert boundary edge limit_prev_pos -> endpoint_prev
          bd.insert(
            iter,
            make_final_or_supported_boundary(iter->limit_prev_anchor, endpoint_prev, prev_blocker));
        }
        // Insert boundary edge endpoint_prev -> V (along obstacle edge)
        bd.insert(iter, make_obstacle_boundary(endpoint_prev, V, *prev_blocker));
      } else {
        // No obstacle edge between V and prev_pos:
        // insert a gap (non-boundary) edge prev_pos -> V with angular sub-range.
        bd.insert(iter,
                  make_portal(iter->prev, V, iter->limit_prev, theta, iter->limit_prev_anchor, V));
        to_minus = true;
      }

      const auto next_blocker = obstacleEdgeBetween(iter->next.position, V.position);
      if (next_blocker) {
        // Obstacle edge between next_pos and V.
        // Compute far-endpoint where ray from source through limit_next_pos hits this edge.
        const auto endpoint_next_position = findIntersection(
          exactPoint(V), exactPoint(iter->next), exact_source, exactPoint(iter->limit_next_anchor));
        if (!endpoint_next_position)
          return false;
        const auto endpoint_next = makeIntersectionEndpoint(*endpoint_next_position, *next_blocker);
        // Insert boundary edge V -> endpoint_next
        auto ni = bd.insert(iter, make_obstacle_boundary(V, endpoint_next, *next_blocker));

        if (exactPoint(iter->limit_next_anchor) != exactPoint(endpoint_next)) {
          // Far-endpoint distinct from limit_next_pos:
          // update iter to boundary edge endpoint_next -> limit_next_pos
          set_final_or_supported_boundary(
            *iter, endpoint_next, iter->limit_next_anchor, next_blocker);
          if (to_minus) {
            iter = std::prev(ni);
          } else {
            iter++;
          }
        } else {
          // Far-endpoint coincides with limit_next_pos: erase iter
          bd.erase(iter++);
          if (to_minus) {
            iter = std::prev(ni);
          }
        }
      } else {
        // No obstacle edge between next_pos and V:
        // update iter to edge V -> next_pos with angular sub-range.
        iter->prev = V;
        iter->limit_prev = theta;
        iter->limit_prev_anchor = V;
        iter->supporting_segment.reset();
        if (to_minus)
          iter--;
      }
    } else if (ray_position == exact_geometry::PortalRayPosition::equal_upper) {
      // CASE 4: theta == limit_next
      // V is collinear with source and next_pos (same angle).
      const auto prev_blocker = obstacleEdgeBetween(V.position, iter->prev.position);
      const auto next_blocker = obstacleEdgeBetween(iter->next.position, V.position);
      if (prev_blocker) {
        // Obstacle edge between V and prev_pos.
        const auto endpoint_prev_position = findIntersection(
          exactPoint(V), exactPoint(iter->prev), exact_source, exactPoint(iter->limit_prev_anchor));
        if (!endpoint_prev_position)
          return false;
        const auto endpoint_prev = makeIntersectionEndpoint(*endpoint_prev_position, *prev_blocker);

        if (exactPoint(endpoint_prev) != exactPoint(iter->limit_prev_anchor)) {
          bd.insert(
            iter,
            make_final_or_supported_boundary(iter->limit_prev_anchor, endpoint_prev, prev_blocker));
        }

        bd.insert(iter, make_obstacle_boundary(endpoint_prev, V, *prev_blocker));

        const auto old_next = iter->next;
        set_final_or_supported_boundary(*iter, V, old_next, next_blocker);
        iter++;
      } else if (next_blocker) {
        // V is collinear with next_pos, but the interior of the angular
        // window remains visible through the edge (prev_pos, V). Continue
        // expanding that portal and preserve the original limit ray so a
        // later constrained edge can produce the far-to-near radial boundary.
        auto ni = bd.insert(iter,
                            make_portal(iter->prev,
                                        V,
                                        iter->limit_prev,
                                        iter->limit_next,
                                        iter->limit_prev_anchor,
                                        iter->limit_next_anchor));
        bd.erase(iter);
        iter = ni;
      } else {
        // No obstacle edge: create gap edge, mark iter as boundary.
        auto ni = bd.insert(
          iter, make_portal(iter->prev, V, iter->limit_prev, theta, iter->limit_prev_anchor, V));

        const auto old_next = iter->next;
        set_final_boundary(*iter, V, old_next);
        // Process the gap edge next (instead of skipping it with iter++).
        // The gap edge ni was inserted before iter; without this, ni would
        // never be expanded, causing visible obstacles behind it to be missed.
        iter = ni;
      }
    } else if (ray_position == exact_geometry::PortalRayPosition::after_upper) {
      // CASE 5: theta > limit_next
      // V is at a larger angle than next_pos.
      const auto blocker = obstacleEdgeBetween(V.position, iter->prev.position);
      if (blocker) {
        // Obstacle edge between V and prev_pos.
        const auto endpoint_prev_position = findIntersection(
          exactPoint(V), exactPoint(iter->prev), exact_source, exactPoint(iter->limit_prev_anchor));
        const auto endpoint_next_position = findIntersection(
          exactPoint(V), exactPoint(iter->prev), exact_source, exactPoint(iter->limit_next_anchor));
        if (!endpoint_prev_position || !endpoint_next_position)
          return false;
        const auto endpoint_prev = makeIntersectionEndpoint(*endpoint_prev_position, *blocker);
        const auto endpoint_next = makeIntersectionEndpoint(*endpoint_next_position, *blocker);

        if (exactPoint(endpoint_prev) != exactPoint(iter->limit_prev_anchor)) {
          bd.insert(
            iter,
            make_final_or_supported_boundary(iter->limit_prev_anchor, endpoint_prev, blocker));
        }

        bd.insert(iter, make_obstacle_boundary(endpoint_prev, endpoint_next, *blocker));

        set_final_or_supported_boundary(*iter, endpoint_next, iter->limit_next_anchor, blocker);
        iter++;
      } else {
        // Symmetric to CASE 1: V lies outside the upper angular limit, so
        // advance the portal endpoint without emitting out-of-window geometry.
        iter->next = V;
      }
    } else {
      return false;
    }
  }

  if (bd.empty())
    return false;
  for (const auto& edge : bd) {
    if (stop_token.poll())
      return false;
    if (!validate_bungiu_edge(edge))
      return false;
    if (exactPoint(edge.prev) == exact_source || exactPoint(edge.next) == exact_source) {
      return false;
    }
  }

  // --- Phase 3: Output ---
  // Sort by exact ray order, while retaining the sweep list position for
  // same-ray visibility discontinuities.  A closed boundary list is cyclic:
  // if one equal-ray group straddles bd.end()/bd.begin(), an ordinary stable
  // sort sees the wrong linear order and can turn far->near into near->far.
  // Repair only that cyclic tie order; all non-equal ray ordering remains the
  // same as the original exact angular sort.
  struct OrderedBoundaryEdge {
    decltype(bd)::iterator edge;
  };
  std::vector<OrderedBoundaryEdge> sorted_edges;
  sorted_edges.reserve(bd.size());
  for (auto it = bd.begin(); it != bd.end(); ++it) {
    if (stop_token.poll())
      return false;
    sorted_edges.push_back(OrderedBoundaryEdge{it});
  }

  if (open_visibility_region) {
    const exact_geometry::Point reference(open_region_start_vertex.first,
                                          open_region_start_vertex.second);
    if (stop_token.poll())
      return false;
    std::stable_sort(sorted_edges.begin(),
                     sorted_edges.end(),
                     [&exact_source, &reference](const auto& first, const auto& second) {
                       return exact_geometry::compareRaysCounterClockwiseFrom(
                                exact_source,
                                reference,
                                exactPoint(first.edge->prev),
                                exactPoint(second.edge->prev)) == CGAL::SMALLER;
                     });
    if (stop_token.poll())
      return false;
  } else {
    const auto compare_closed_rays = [&exact_source](const auto& first, const auto& second) {
      return exact_geometry::compareRaysLikeAtan2(
        exact_source, exactPoint(first.edge->prev), exactPoint(second.edge->prev));
    };
    if (stop_token.poll())
      return false;
    std::stable_sort(sorted_edges.begin(),
                     sorted_edges.end(),
                     [&compare_closed_rays](const auto& first, const auto& second) {
                       return compare_closed_rays(first, second) == CGAL::SMALLER;
                     });
    if (stop_token.poll())
      return false;

    size_t group_begin = 0;
    while (group_begin < sorted_edges.size()) {
      if (stop_token.poll())
        return false;
      size_t group_end = group_begin + 1;
      while (group_end < sorted_edges.size() &&
             compare_closed_rays(sorted_edges[group_begin], sorted_edges[group_end]) ==
               CGAL::EQUAL) {
        if (stop_token.poll())
          return false;
        ++group_end;
      }

      if (group_end - group_begin > 1) {
        const size_t group_size = group_end - group_begin;
        const auto rotation_satisfies_radial_edges = [&](size_t rotation) {
          if (stop_token.poll())
            return false;
          const auto rotated_position = [&](size_t local_index) {
            return (local_index + group_size - rotation) % group_size;
          };

          for (size_t local_from = 0; local_from < group_size; ++local_from) {
            if (stop_token.poll())
              return false;
            const auto& from = sorted_edges[group_begin + local_from];
            const exact_geometry::Point from_point = exactPoint(from.edge->prev);
            const exact_geometry::Point next_point = exactPoint(from.edge->next);
            if (from_point == next_point ||
                !exact_geometry::isSameDirectedRay(exact_source, from_point, next_point)) {
              continue;
            }

            bool has_group_target = false;
            bool direction_satisfied = false;
            for (size_t local_to = 0; local_to < group_size; ++local_to) {
              if (stop_token.poll())
                return false;
              const auto& target = sorted_edges[group_begin + local_to];
              if (exactPoint(target.edge->prev) != next_point ||
                  !(target.edge->prev.support == from.edge->next.support)) {
                continue;
              }
              has_group_target = true;
              if (rotated_position(local_from) < rotated_position(local_to)) {
                direction_satisfied = true;
                break;
              }
            }
            if (!has_group_target || !direction_satisfied)
              return false;
          }
          return true;
        };

        // Keep the original stable order whenever it respects every explicit
        // same-ray boundary edge.  Otherwise accept only one unambiguous
        // cyclic rotation; guessing between several candidates could silently
        // reconnect the boundary on the wrong side of an occluder.
        if (!rotation_satisfies_radial_edges(0)) {
          size_t valid_rotation = 0;
          size_t valid_rotation_count = 0;
          for (size_t rotation = 1; rotation < group_size; ++rotation) {
            if (stop_token.poll())
              return false;
            if (rotation_satisfies_radial_edges(rotation)) {
              valid_rotation = rotation;
              ++valid_rotation_count;
            }
          }
          if (valid_rotation_count != 1)
            return false;
          std::rotate(
            sorted_edges.begin() + static_cast<std::ptrdiff_t>(group_begin),
            sorted_edges.begin() + static_cast<std::ptrdiff_t>(group_begin + valid_rotation),
            sorted_edges.begin() + static_cast<std::ptrdiff_t>(group_end));
        }
      }
      group_begin = group_end;
    }
  }

  for (const auto& ordered : sorted_edges) {
    if (stop_token.poll()) {
      visibility_region.clear();
      return false;
    }
    const auto it = ordered.edge;
    if (!std::isfinite(it->prev.position.first) || !std::isfinite(it->prev.position.second)) {
      visibility_region.clear();
      return false;
    }
    visibility_region.emplace_back(it->prev);
  }

  if (open_visibility_region)
    visibility_region.emplace_back(makeEndpoint(open_region_closing_vertex));

  if (visibility_region.size() < 2) {
    visibility_region.clear();
    return false;
  }

  std::string validation_error;
  if (!validateVisibilityRegionImpl(visibility_region, stop_token, &validation_error)) {
    visibility_region.clear();
    return false;
  }

  std::string context_error;
  const auto context = makeVisibilityGeometryContext(
    *this, source, exact_source, obstacle_vertex_source, context_error);
  if (!context ||
      detail::validateVisibilityGeometry(
        visibility_region, *context, stop_token, &validation_error) != OperationStatus::success) {
    visibility_region.clear();
    return false;
  }

  return true;
}

bool Polymap::calculateVisibilityRegion(int x,
                                        int y,
                                        std::vector<std::pair<double, double>>& result_V,
                                        std::vector<std::pair<int, int>>& topo_V) {
  result_V.clear();
  topo_V.clear();

  VisibilityRegion visibility_region;
  if (!calculateVisibilityRegion(x, y, visibility_region))
    return false;

  result_V.reserve(visibility_region.size());
  topo_V.reserve(visibility_region.size());
  for (const auto& endpoint : visibility_region) {
    result_V.emplace_back(endpoint.position);
    if (const auto* vertex = std::get_if<ObstacleVertexId>(&endpoint.support))
      topo_V.emplace_back(toLegacyTopology(*vertex));
    else
      topo_V.emplace_back(-1, -1);
  }
  return true;
}

BoundaryEndpoint Polymap::makeEndpoint(const Point2d& position) const {
  BoundaryEndpoint endpoint;
  endpoint.position = position;
  endpoint.exact_position = exact_geometry::Point(position.first, position.second);
  const auto topology = locateVertex(position.first, position.second);
  if (isValidTopology(topology)) {
    const auto& vertex = obs_[static_cast<size_t>(topology.first)]
                           .ordered_vertices_[static_cast<size_t>(topology.second)];
    endpoint.position = {static_cast<double>(vertex.first), static_cast<double>(vertex.second)};
    endpoint.exact_position = exact_geometry::Point(vertex.first, vertex.second);
    endpoint.support = fromLegacyTopology(topology);
  }
  return endpoint;
}

BoundaryEndpoint Polymap::makeEndpoint(const std::pair<int, int>& position) const {
  BoundaryEndpoint endpoint;
  endpoint.position = {static_cast<double>(position.first), static_cast<double>(position.second)};
  endpoint.exact_position = exact_geometry::Point(position.first, position.second);
  const auto topology = locateVertex(position.first, position.second);
  if (isValidTopology(topology))
    endpoint.support = fromLegacyTopology(topology);
  return endpoint;
}

BoundaryEndpoint Polymap::makeIntersectionEndpoint(
  const exact_geometry::Point& position, const DirectedObstacleEdge& supporting_edge) const {
  BoundaryEndpoint endpoint;
  endpoint.position = projectExactPoint(position);
  endpoint.exact_position = position;
  endpoint.support = supporting_edge;

  const auto from_topology = toLegacyTopology(supporting_edge.from);
  const auto to_topology = toLegacyTopology(supporting_edge.to);
  if (!isValidTopology(from_topology) || !isValidTopology(to_topology) ||
      !areConsecutive(from_topology, to_topology)) {
    return endpoint;
  }

  const auto& from = obs_[static_cast<size_t>(supporting_edge.from.obstacle)]
                       .ordered_vertices_[static_cast<size_t>(supporting_edge.from.vertex)];
  const auto& to = obs_[static_cast<size_t>(supporting_edge.to.obstacle)]
                     .ordered_vertices_[static_cast<size_t>(supporting_edge.to.vertex)];
  const exact_geometry::Point exact_from(from.first, from.second);
  const exact_geometry::Point exact_to(to.first, to.second);
  if (position == exact_from) {
    endpoint.position = {static_cast<double>(from.first), static_cast<double>(from.second)};
    endpoint.exact_position = exact_from;
    endpoint.support = supporting_edge.from;
  } else if (position == exact_to) {
    endpoint.position = {static_cast<double>(to.first), static_cast<double>(to.second)};
    endpoint.exact_position = exact_to;
    endpoint.support = supporting_edge.to;
  } else {
    endpoint.support = supporting_edge;
  }
  return endpoint;
}

std::optional<DirectedObstacleEdge> Polymap::obstacleEdgeBetween(const Point2d& from,
                                                                 const Point2d& to) const {
  const auto from_topology = locateVertex(from.first, from.second);
  const auto to_topology = locateVertex(to.first, to.second);
  if (!areConsecutive(from_topology, to_topology))
    return std::nullopt;
  return DirectedObstacleEdge{fromLegacyTopology(from_topology), fromLegacyTopology(to_topology)};
}

bool Polymap::isAnObstacleEdge(const Point2d& from, const Point2d& to) const {
  return obstacleEdgeBetween(from, to).has_value();
}

bool Polymap::pointLiesOnObstacleEdge(const BoundaryEndpoint& point,
                                      const DirectedObstacleEdge& edge) const {
  if (!std::isfinite(point.position.first) || !std::isfinite(point.position.second))
    return false;

  const auto from_topology = toLegacyTopology(edge.from);
  const auto to_topology = toLegacyTopology(edge.to);
  if (!areConsecutive(from_topology, to_topology))
    return false;

  const auto& from = obs_[static_cast<size_t>(edge.from.obstacle)]
                       .ordered_vertices_[static_cast<size_t>(edge.from.vertex)];
  const auto& to = obs_[static_cast<size_t>(edge.to.obstacle)]
                     .ordered_vertices_[static_cast<size_t>(edge.to.vertex)];

  const exact_geometry::Point exact_from(from.first, from.second);
  const exact_geometry::Point exact_to(to.first, to.second);
  const exact_geometry::Point exact_point = exactPoint(point);
  return CGAL::orientation(exact_from, exact_point, exact_to) == CGAL::COLLINEAR &&
         CGAL::collinear_are_ordered_along_line(exact_from, exact_point, exact_to);
}

bool Polymap::validateBoundaryEndpoint(const BoundaryEndpoint& endpoint, std::string* error) const {
  const StopToken no_stop;
  return validateBoundaryEndpointImpl(endpoint, no_stop, error);
}

OperationStatus Polymap::validateBoundaryEndpoint(const BoundaryEndpoint& endpoint,
                                                  const StopToken& stop_token,
                                                  std::string* error) const {
  if (error)
    error->clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (validateBoundaryEndpointImpl(endpoint, stop_token, error)) {
    if (stop_token.poll()) {
      if (error)
        error->clear();
      return OperationStatus::stopped;
    }
    return OperationStatus::success;
  }
  if (stop_token.poll()) {
    if (error)
      error->clear();
    return OperationStatus::stopped;
  }
  return OperationStatus::failure;
}

bool Polymap::validateBoundaryEndpointImpl(const BoundaryEndpoint& endpoint,
                                           const StopToken& stop_token,
                                           std::string* error) const {
  if (error)
    error->clear();
  const auto fail = [&](const std::string& message) {
    if (error)
      *error = message;
    return false;
  };

  if (stop_token.poll())
    return false;

  if (!std::isfinite(endpoint.position.first) || !std::isfinite(endpoint.position.second)) {
    return fail("Boundary endpoint position is not finite");
  }
  if (endpoint.exact_position && endpoint.position != projectExactPoint(*endpoint.exact_position)) {
    return fail("Boundary endpoint double projection does not match its exact position");
  }

  if (std::holds_alternative<std::monostate>(endpoint.support)) {
    if (isValidTopology(locateVertex(endpoint.position.first, endpoint.position.second)))
      return fail("Obstacle vertex endpoint has no exact-vertex support");
    for (size_t obstacle_index = 0; obstacle_index < obs_.size(); ++obstacle_index) {
      if (stop_token.poll())
        return false;
      const auto& vertices = obs_[obstacle_index].ordered_vertices_;
      for (size_t vertex_index = 0; vertex_index < vertices.size(); ++vertex_index) {
        if (stop_token.poll())
          return false;
        const size_t next_index = (vertex_index + 1) % vertices.size();
        const DirectedObstacleEdge edge{
          ObstacleVertexId{static_cast<int>(obstacle_index), static_cast<int>(vertex_index)},
          ObstacleVertexId{static_cast<int>(obstacle_index), static_cast<int>(next_index)}};
        if (pointLiesOnObstacleEdge(endpoint, edge))
          return fail("Obstacle-edge endpoint has no supporting-edge metadata");
      }
    }
    return true;
  }

  if (const auto* vertex = std::get_if<ObstacleVertexId>(&endpoint.support)) {
    const auto topology = toLegacyTopology(*vertex);
    if (!isValidTopology(topology))
      return fail("Boundary endpoint references an invalid obstacle vertex");
    const auto& expected = obs_[static_cast<size_t>(vertex->obstacle)]
                             .ordered_vertices_[static_cast<size_t>(vertex->vertex)];
    constexpr double coordinate_tolerance = 1e-9;
    if (std::abs(endpoint.position.first - static_cast<double>(expected.first)) >
          coordinate_tolerance ||
        std::abs(endpoint.position.second - static_cast<double>(expected.second)) >
          coordinate_tolerance) {
      return fail("Boundary endpoint position does not match its obstacle vertex");
    }
    if (exactPoint(endpoint) != exact_geometry::Point(expected.first, expected.second))
      return fail("Boundary endpoint exact position does not match its obstacle vertex");
    return true;
  }

  const auto* edge = std::get_if<DirectedObstacleEdge>(&endpoint.support);
  if (!edge || !pointLiesOnObstacleEdge(endpoint, *edge))
    return fail("Boundary endpoint is not on its supporting obstacle edge");

  const auto& from = obs_[static_cast<size_t>(edge->from.obstacle)]
                       .ordered_vertices_[static_cast<size_t>(edge->from.vertex)];
  const auto& to = obs_[static_cast<size_t>(edge->to.obstacle)]
                     .ordered_vertices_[static_cast<size_t>(edge->to.vertex)];
  const bool is_from = exactPoint(endpoint) == exact_geometry::Point(from.first, from.second);
  const bool is_to = exactPoint(endpoint) == exact_geometry::Point(to.first, to.second);
  if (is_from || is_to)
    return fail("Obstacle-edge endpoint support was not canonicalized to an exact vertex");
  return true;
}

bool Polymap::validateVisibilityRegion(const VisibilityRegion& visibility_region,
                                       std::string* error) const {
  const StopToken no_stop;
  return validateVisibilityRegionImpl(visibility_region, no_stop, error);
}

OperationStatus Polymap::validateVisibilityRegion(const VisibilityRegion& visibility_region,
                                                  const StopToken& stop_token,
                                                  std::string* error) const {
  if (error)
    error->clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (validateVisibilityRegionImpl(visibility_region, stop_token, error)) {
    if (stop_token.poll()) {
      if (error)
        error->clear();
      return OperationStatus::stopped;
    }
    return OperationStatus::success;
  }
  if (stop_token.poll()) {
    if (error)
      error->clear();
    return OperationStatus::stopped;
  }
  return OperationStatus::failure;
}

bool Polymap::validateVisibilityRegionImpl(const VisibilityRegion& visibility_region,
                                           const StopToken& stop_token,
                                           std::string* error) const {
  if (error)
    error->clear();
  if (stop_token.poll())
    return false;
  if (visibility_region.size() < 2) {
    if (error)
      *error = "Visibility region has fewer than two vertices";
    return false;
  }

  for (size_t i = 0; i < visibility_region.size(); ++i) {
    if (stop_token.poll())
      return false;
    std::string endpoint_error;
    if (!validateBoundaryEndpointImpl(visibility_region[i], stop_token, &endpoint_error)) {
      if (error) {
        *error = "Visibility endpoint " + std::to_string(i) + " is invalid: " + endpoint_error;
      }
      return false;
    }
  }
  return true;
}

bool Polymap::boundarySupportsConsecutive(const BoundaryEndpoint& prev,
                                          const BoundaryEndpoint& next) const {
  const StopToken no_stop;
  bool supports = false;
  if (!boundarySupportsConsecutiveImpl(prev, next, supports, no_stop))
    return false;
  return supports;
}

OperationStatus Polymap::boundarySupportsConsecutive(const BoundaryEndpoint& prev,
                                                     const BoundaryEndpoint& next,
                                                     bool& supports,
                                                     const StopToken& stop_token) const {
  supports = false;
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (boundarySupportsConsecutiveImpl(prev, next, supports, stop_token)) {
    if (stop_token.poll()) {
      supports = false;
      return OperationStatus::stopped;
    }
    return OperationStatus::success;
  }
  supports = false;
  return stop_token.poll() ? OperationStatus::stopped : OperationStatus::failure;
}

bool Polymap::boundarySupportsConsecutiveImpl(const BoundaryEndpoint& prev,
                                              const BoundaryEndpoint& next,
                                              bool& supports,
                                              const StopToken& stop_token) const {
  supports = false;
  if (!validateBoundaryEndpointImpl(prev, stop_token, nullptr) ||
      !validateBoundaryEndpointImpl(next, stop_token, nullptr))
    return false;

  const auto* prev_vertex = std::get_if<ObstacleVertexId>(&prev.support);
  const auto* next_vertex = std::get_if<ObstacleVertexId>(&next.support);
  const auto* prev_edge = std::get_if<DirectedObstacleEdge>(&prev.support);
  const auto* next_edge = std::get_if<DirectedObstacleEdge>(&next.support);

  if (prev_vertex && next_vertex) {
    supports = areConsecutive(toLegacyTopology(*prev_vertex), toLegacyTopology(*next_vertex));
    return true;
  }
  if (prev_vertex && next_edge) {
    supports = *prev_vertex == next_edge->from;
    return true;
  }
  if (prev_edge && next_vertex) {
    supports = prev_edge->to == *next_vertex;
    return true;
  }
  if (prev_edge && next_edge) {
    if (!(*prev_edge == *next_edge)) {
      supports = prev_edge->to == next_edge->from;
      return true;
    }

    const auto from_topology = toLegacyTopology(prev_edge->from);
    const auto to_topology = toLegacyTopology(prev_edge->to);
    if (!areConsecutive(from_topology, to_topology))
      return true;
    const auto& from = obs_[static_cast<size_t>(prev_edge->from.obstacle)]
                         .ordered_vertices_[static_cast<size_t>(prev_edge->from.vertex)];
    const auto& to = obs_[static_cast<size_t>(prev_edge->to.obstacle)]
                       .ordered_vertices_[static_cast<size_t>(prev_edge->to.vertex)];
    const exact_geometry::Point exact_from(from.first, from.second);
    const exact_geometry::Point exact_to(to.first, to.second);
    const auto direction = exact_to - exact_from;
    const auto length_squared = direction.squared_length();
    if (length_squared == exact_geometry::FT(0))
      return true;
    const auto parameter_numerator = [&](const BoundaryEndpoint& endpoint) {
      return CGAL::scalar_product(exactPoint(endpoint) - exact_from, direction);
    };
    supports = parameter_numerator(prev) <= parameter_numerator(next);
    return true;
  }
  return true;
}

bool Polymap::isInTri(int x1, int y1, int x2, int y2, int x3, int y3, double x, double y) {
  const exact_geometry::Kernel::Triangle_2 triangle(
    exact_geometry::Point(x1, y1), exact_geometry::Point(x2, y2), exact_geometry::Point(x3, y3));
  return triangle.bounded_side(exact_geometry::Point(x, y)) != CGAL::ON_UNBOUNDED_SIDE;
}

void Polymap::clearCGALRelatedState() {
  cdt_ready_ = false;
  cdt_.clear();
  facets_.clear();
  triangle_faces_.clear();
  triangle_edges_.clear();
  cdt_table_.clear();
  cdt_ver_num_ = 0;
  visibility_storage_.clear();
}

void Polymap::clearStoppedConstructionState() {
  solution_exist_ = false;
  no_path_ = false;
  construction_stopped_ = true;
  construction_error_.clear();
  obs_.clear();
  std::fill(vertices_location_x_flat_.begin(), vertices_location_x_flat_.end(), -1);
  std::fill(vertices_location_y_flat_.begin(), vertices_location_y_flat_.end(), -1);
  clearCGALRelatedState();
}

bool Polymap::constructCGALRelated(CdtValidator validator, std::string& error) {
  const StopToken no_stop;
  return constructCGALRelatedImpl(validator, error, no_stop);
}

OperationStatus Polymap::constructCGALRelated(CdtValidator validator,
                                              std::string& error,
                                              const StopToken& stop_token) {
  error.clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (constructCGALRelatedImpl(validator, error, stop_token))
    return OperationStatus::success;
  if (stop_token.poll()) {
    error.clear();
    return OperationStatus::stopped;
  }
  return OperationStatus::failure;
}

bool Polymap::constructCGALRelatedImpl(CdtValidator validator,
                                       std::string& error,
                                       const StopToken& stop_token) {
  using CDT = constrained_delaunay_triangulation::CDT;

  error.clear();
  clearCGALRelatedState();

  const auto fail = [&](std::string message) {
    clearCGALRelatedState();
    error = std::move(message);
    return false;
  };

  if (validator == nullptr)
    return fail("CDT validation function is unavailable");
  if (stop_token.poll())
    return false;

  try {
    CDT candidate_cdt;
    for (const auto& ob : obs_) {
      if (stop_token.poll())
        return false;
      for (auto it = ob.ordered_vertices_.begin(); it != ob.ordered_vertices_.end(); ++it) {
        if (stop_token.poll())
          return false;
        auto next = std::next(it);
        if (next == ob.ordered_vertices_.end())
          next = ob.ordered_vertices_.begin();
        candidate_cdt.insert_constraint(
          constrained_delaunay_triangulation::Point(it->first, it->second),
          constrained_delaunay_triangulation::Point(next->first, next->second));
        if (stop_token.poll())
          return false;
      }
    }

    if (stop_token.poll())
      return false;
    const bool valid_cdt = validator(candidate_cdt);
    if (stop_token.poll())
      return false;
    if (!valid_cdt)
      return fail("CGAL reported that the constrained triangulation is invalid");

    std::vector<std::vector<std::pair<int, int>>> candidate_facets;
    std::unordered_map<long long, int> candidate_table;
    const int candidate_vertex_count = static_cast<int>(candidate_cdt.number_of_vertices());
    const long long map_area = static_cast<long long>(xsize_) * static_cast<long long>(ysize_);
    const auto vertex_key = [&](const std::pair<int, int>& vertex) {
      return static_cast<long long>(vertex.first) +
             static_cast<long long>(vertex.second) * static_cast<long long>(xsize_);
    };
    const auto directed_edge_key = [&](const std::pair<int, int>& from,
                                       const std::pair<int, int>& to) {
      return vertex_key(from) + vertex_key(to) * map_area;
    };

    int count = 0;
    for (auto fit = candidate_cdt.finite_faces_begin(); fit != candidate_cdt.finite_faces_end();
         ++fit) {
      if (stop_token.poll())
        return false;
      candidate_facets.emplace_back(
        std::vector<std::pair<int, int>>{{static_cast<int>(fit->vertex(0)->point().x()),
                                          static_cast<int>(fit->vertex(0)->point().y())},
                                         {static_cast<int>(fit->vertex(1)->point().x()),
                                          static_cast<int>(fit->vertex(1)->point().y())},
                                         {static_cast<int>(fit->vertex(2)->point().x()),
                                          static_cast<int>(fit->vertex(2)->point().y())}});
      const auto& facet = candidate_facets.back();
      fit->info().stable_id = count;
      fit->info().is_free = false;
      candidate_table[directed_edge_key(facet[0], facet[1])] = count;
      candidate_table[directed_edge_key(facet[1], facet[2])] = count;
      candidate_table[directed_edge_key(facet[2], facet[0])] = count;
      ++count;
    }

    if (stop_token.poll())
      return false;
    cdt_ = std::move(candidate_cdt);
    facets_ = std::move(candidate_facets);
    cdt_table_ = std::move(candidate_table);
    cdt_ver_num_ = candidate_vertex_count;
    cdt_ready_ = true;
    const OperationStatus triangle_environment_status =
      buildTriangleEnvironment(error, stop_token);
    if (triangle_environment_status == OperationStatus::stopped)
      return false;
    if (triangle_environment_status == OperationStatus::failure)
      return fail(error.empty() ? "Could not build the free-triangle environment" : error);
    return true;
  } catch (const CDT::Intersection_of_constraints_exception&) {
    return fail("Obstacle constraints intersect, overlap, or form an unsupported T-junction");
  } catch (const CGAL::Failure_exception& exception) {
    return fail("CGAL rejected the constrained triangulation: " + std::string(exception.what()));
  } catch (const std::exception& exception) {
    return fail("Unexpected error while constructing the CDT: " + std::string(exception.what()));
  } catch (...) {
    return fail("Unknown error while constructing the CDT");
  }
}

bool Polymap::getPolyObstacles(int start_x, int start_y, int goal_x, int goal_y) {
  const StopToken no_stop;
  return getPolyObstacles(start_x, start_y, goal_x, goal_y, no_stop) == OperationStatus::success;
}

OperationStatus Polymap::getPolyObstacles(
  int start_x, int start_y, int goal_x, int goal_y, const StopToken& stop_token) {
  return getPolyObstacles(
    start_x,
    start_y,
    std::vector<PolymapEndpoint>{
      {goal_x, goal_y, {static_cast<double>(goal_x), static_cast<double>(goal_y)}}},
    stop_token);
}

OperationStatus Polymap::getPolyObstacles(int start_x,
                                          int start_y,
                                          const std::vector<PolymapEndpoint>& goals,
                                          const StopToken& stop_token) {
  return getPolyObstacles(start_x, start_y, goals, stop_token, std::nullopt);
}

OperationStatus Polymap::getPolyObstacles(
  int start_x,
  int start_y,
  const std::vector<PolymapEndpoint>& goals,
  const StopToken& stop_token,
  std::optional<size_t> max_raw_contour_vertices) {
  if (stop_token.poll())
    return OperationStatus::stopped;
  construction_error_.clear();
  if (getPolyObstaclesImpl(
        start_x, start_y, goals, stop_token, max_raw_contour_vertices)) {
    // obs_ now contains a newly extracted, unsimplified contour set.  Any
    // topology registry, CDT, facet table or visibility cache from a previous
    // build refers to different vertex indices and must not remain usable.
    solution_exist_ = true;
    construction_stopped_ = false;
    construction_error_.clear();
    std::fill(vertices_location_x_flat_.begin(), vertices_location_x_flat_.end(), -1);
    std::fill(vertices_location_y_flat_.begin(), vertices_location_y_flat_.end(), -1);
    clearCGALRelatedState();
    return OperationStatus::success;
  }
  if (stop_token.poll())
    return OperationStatus::stopped;
  // A normal refresh failure invalidates every object derived from the old
  // obstacle set; keeping a ready CDT/cache next to an empty obs_ vector would
  // expose stale geometry.  A cooperative stop is different: the previously
  // committed state remains untouched.
  solution_exist_ = false;
  construction_stopped_ = false;
  obs_.clear();
  std::fill(vertices_location_x_flat_.begin(), vertices_location_x_flat_.end(), -1);
  std::fill(vertices_location_y_flat_.begin(), vertices_location_y_flat_.end(), -1);
  clearCGALRelatedState();
  return OperationStatus::failure;
}

bool Polymap::getPolyObstaclesImpl(
  int start_x,
  int start_y,
  int goal_x,
  int goal_y,
  const StopToken& stop_token,
  std::optional<size_t> max_raw_contour_vertices) {
  return getPolyObstaclesImpl(
    start_x,
    start_y,
    std::vector<PolymapEndpoint>{
      {goal_x, goal_y, {static_cast<double>(goal_x), static_cast<double>(goal_y)}}},
    stop_token,
    max_raw_contour_vertices);
}

bool Polymap::getPolyObstaclesImpl(int start_x,
                                   int start_y,
                                   const std::vector<PolymapEndpoint>& goals,
                                   const StopToken& stop_token,
                                   std::optional<size_t> max_raw_contour_vertices) {
  if (stop_token.poll())
    return false;
  if (xsize_ <= 0 || ysize_ <= 0)
    return false;
  const size_t width = static_cast<size_t>(xsize_);
  const size_t height = static_cast<size_t>(ysize_);
  if (width > std::numeric_limits<size_t>::max() / height)
    return false;
  const size_t cell_count = width * height;
  if (cell_count > static_cast<size_t>(std::numeric_limits<int>::max()) ||
      data_.size() != cell_count) {
    return false;
  }
  const auto in_bounds = [this](int x, int y) {
    return x >= 0 && y >= 0 && x < xsize_ && y < ysize_;
  };
  if (!in_bounds(start_x, start_y) || goals.empty())
    return false;
  for (const auto& goal : goals) {
    if (!in_bounds(goal.cell_x, goal.cell_y))
      return false;
  }

  std::vector<Obs> candidate_obstacles;
  unsigned int nx = xsize_;
  unsigned int ny = ysize_;

  std::vector<int> mask(nx * ny, 0);

  std::unordered_map<UndirectedGridEdgeKey, DirectedGridEdge, UndirectedGridEdgeKeyHash> edges;
  const auto add_directed_edge = [&](int from, int to) {
    const auto insertion =
      edges.emplace(makeUndirectedGridEdgeKey(from, to), DirectedGridEdge{from, to});
    if (max_raw_contour_vertices && edges.size() > *max_raw_contour_vertices) {
      if (insertion.second)
        edges.erase(insertion.first);
      construction_error_ =
        "Reference-shortening raw contour exceeds the vertex budget " +
        std::to_string(*max_raw_contour_vertices);
      return false;
    }
    return true;
  };
  std::stack<int> Q;
  Q.push(start_x + start_y * nx);
  while (!Q.empty()) {
    if (stop_token.poll())
      return false;
    int cur = Q.top();
    Q.pop();
    int x = cur % nx;
    int y = (cur - x) / nx;
    if (data_[cur] != 0 || mask[x + y * nx] != 0)
      continue;

    if (cur - 1 >= 0 && data_[cur - 1] != 0 && !add_directed_edge(cur, cur + nx))
      return false;
    if (cur + 1 < static_cast<int>(nx * ny) && data_[cur + 1] != 0 &&
        !add_directed_edge(cur + nx + 1, cur + 1)) {
      return false;
    }
    if (cur - static_cast<int>(nx) >= 0 && data_[cur - nx] != 0 &&
        !add_directed_edge(cur + 1, cur)) {
      return false;
    }
    if (cur + static_cast<int>(nx) < static_cast<int>(nx * ny) && data_[cur + nx] != 0 &&
        !add_directed_edge(cur + nx, cur + nx + 1)) {
      return false;
    }

    mask[x + y * nx] = 1;
    if (x > 0 && data_[cur - 1] == 0)
      Q.push(cur - 1);
    if (x < static_cast<int>(nx) - 1 && data_[cur + 1] == 0)
      Q.push(cur + 1);
    if (y > 0 && data_[cur - nx] == 0)
      Q.push(cur - nx);
    if (y < static_cast<int>(ny) - 1 && data_[cur + nx] == 0)
      Q.push(cur + nx);
  }

  for (const auto& goal : goals) {
    if (mask[static_cast<size_t>(goal.cell_x) + static_cast<size_t>(goal.cell_y) * nx] == 0)
      return false;
  }

  std::list<std::pair<int, int>> boundary_points;
  while (!edges.empty()) {
    if (stop_token.poll())
      return false;
    candidate_obstacles.emplace_back(Obs());
    boundary_points.clear();

    auto first_iter = edges.begin();
    const int prev = first_iter->second.from;
    const int next = first_iter->second.to;
    int prev_x = prev % nx;
    int prev_y = (prev - prev_x) / nx;
    int next_x = next % nx;
    int next_y = (next - next_x) / nx;

    boundary_points.emplace_back(prev_x, prev_y);
    boundary_points.emplace_back(next_x, next_y);

    auto biter = boundary_points.end();
    int x, y;
    while (1) {
      if (stop_token.poll())
        return false;
      x = boundary_points.back().first;
      y = boundary_points.back().second;
      int cur = x + y * nx;

      int lb_free = 0, lt_free = 0, rb_free = 0, rt_free = 0;
      if (x > 0 && y > 0 && data_[cur - nx - 1] == 0)
        lb_free = 1;
      if (x > 0 && data_[cur - 1] == 0)
        lt_free = 1;
      if (y > 0 && data_[cur - nx] == 0)
        rb_free = 1;
      if (data_[cur] == 0)
        rt_free = 1;

      int num = lb_free * 8 + lt_free * 4 + rb_free * 2 + rt_free;

      switch (num) {
      case 1:
        boundary_points.emplace_back(x, y + 1);
        break;
      case 2:
        boundary_points.emplace_back(x + 1, y);
        break;
      case 3:
        boundary_points.emplace_back(x, y + 1);
        break;
      case 4:
        boundary_points.emplace_back(x - 1, y);
        break;
      case 5:
        boundary_points.emplace_back(x - 1, y);
        break;
      case 6:
        biter = std::prev(boundary_points.end(), 2);
        if ((*biter).first == x && (*biter).second == y - 1)
          boundary_points.emplace_back(x + 1, y);
        else if ((*biter).first == x && (*biter).second == y + 1)
          boundary_points.emplace_back(x - 1, y);
        break;
      case 7:
        boundary_points.emplace_back(x - 1, y);
        break;
      case 8:
        boundary_points.emplace_back(x, y - 1);
        break;
      case 9:
        biter = std::prev(boundary_points.end(), 2);
        if ((*biter).first == x + 1 && (*biter).second == y)
          boundary_points.emplace_back(x, y + 1);
        else if ((*biter).first == x - 1 && (*biter).second == y)
          boundary_points.emplace_back(x, y - 1);
        break;
      case 10:
        boundary_points.emplace_back(x + 1, y);
        break;
      case 11:
        boundary_points.emplace_back(x, y + 1);
        break;
      case 12:
        boundary_points.emplace_back(x, y - 1);
        break;
      case 13:
        boundary_points.emplace_back(x, y - 1);
        break;
      case 14:
        boundary_points.emplace_back(x + 1, y);
        break;
      }

      if (boundary_points.back().first == boundary_points.front().first &&
          boundary_points.back().second == boundary_points.front().second) {
        boundary_points.pop_back();
        break;
      }
    }

    for (auto& bp : boundary_points) {
      if (stop_token.poll())
        return false;
      candidate_obstacles.back().ordered_vertices_.emplace_back(bp);
    }

    for (auto it = candidate_obstacles.back().ordered_vertices_.begin();
         it != candidate_obstacles.back().ordered_vertices_.end();
         ++it) {
      if (stop_token.poll())
        return false;
      int curr = it->first + it->second * nx;
      int nxt;
      if (it == candidate_obstacles.back().ordered_vertices_.end() - 1)
        nxt = candidate_obstacles.back().ordered_vertices_.front().first +
              candidate_obstacles.back().ordered_vertices_.front().second * nx;
      else
        nxt = std::next(it)->first + std::next(it)->second * nx;
      edges.erase(makeUndirectedGridEdgeKey(curr, nxt));
    }
  }
  if (stop_token.poll())
    return false;
  obs_ = std::move(candidate_obstacles);
  return true;
}

bool Polymap::getVisibilityRegion(int start_x, int start_y, VisibilityRegion& visibility_region) {
  return getVisibilityRegion(start_x, start_y, visibility_region, nullptr);
}

bool Polymap::getVisibilityRegion(int start_x,
                                  int start_y,
                                  VisibilityRegion& visibility_region,
                                  std::string* error) {
  const StopToken no_stop;
  return getVisibilityRegionImpl(start_x, start_y, visibility_region, no_stop, error);
}

OperationStatus Polymap::getVisibilityRegion(int start_x,
                                             int start_y,
                                             VisibilityRegion& visibility_region,
                                             const StopToken& stop_token,
                                             std::string* error) {
  visibility_region.clear();
  if (error)
    error->clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (getVisibilityRegionImpl(start_x, start_y, visibility_region, stop_token, error)) {
    return OperationStatus::success;
  }
  visibility_region.clear();
  if (stop_token.poll()) {
    if (error)
      error->clear();
    return OperationStatus::stopped;
  }
  return OperationStatus::failure;
}

bool Polymap::getVisibilityRegionImpl(int start_x,
                                      int start_y,
                                      VisibilityRegion& visibility_region,
                                      const StopToken& stop_token,
                                      std::string* error) {
  visibility_region.clear();
  if (error)
    error->clear();

  const auto fail = [&](const std::string& message) {
    visibility_region.clear();
    if (error)
      *error = message;
    return false;
  };

  if (stop_token.poll())
    return false;

  if (!solution_exist_)
    return fail("Visibility region is unavailable because the map has no solution");
  if (!cdt_ready_) {
    std::string message =
      "Visibility region is unavailable because map geometry construction failed";
    if (!construction_error_.empty())
      message += ": " + construction_error_;
    return fail(message);
  }
  if (start_x < 0 || start_y < 0 || start_x >= xsize_ || start_y >= ysize_)
    return fail("Visibility source is outside map bounds");

  const auto validate_for_source = [&](const VisibilityRegion& region,
                                       std::string& validation_error) {
    if (!validateVisibilityRegionImpl(region, stop_token, &validation_error))
      return false;

    std::string context_error;
    const auto context = makeVisibilityGeometryContext(*this, start_x, start_y, context_error);
    if (!context) {
      validation_error = std::move(context_error);
      return false;
    }
    return detail::validateVisibilityGeometry(region, *context, stop_token, &validation_error) ==
           OperationStatus::success;
  };

  const int index = start_x + start_y * xsize_;
  auto cached = visibility_storage_.find(index);
  if (cached != visibility_storage_.end()) {
    std::string validation_error;
    if (!validate_for_source(cached->second.vertices, validation_error)) {
      if (stop_token.poll())
        return false;
      visibility_storage_.erase(cached);
      return fail("Cached " + validation_error);
    }
    if (stop_token.poll())
      return false;
    VisibilityRegion candidate_output = cached->second.vertices;
    if (stop_token.poll())
      return false;
    visibility_region = std::move(candidate_output);
    return true;
  }

  VisibilityCacheEntry entry;
  if (!calculateVisibilityRegionImpl(start_x, start_y, entry.vertices, stop_token))
    return fail("Visibility region calculation failed");

  std::string validation_error;
  if (!validate_for_source(entry.vertices, validation_error))
    return fail(validation_error);

  if (stop_token.poll())
    return false;

  VisibilityRegion candidate_output = entry.vertices;
  if (stop_token.poll())
    return false;

  const auto insertion = visibility_storage_.emplace(index, std::move(entry));
  if (stop_token.poll()) {
    if (insertion.second)
      visibility_storage_.erase(insertion.first);
    return false;
  }
  visibility_region = std::move(candidate_output);
  return true;
}

bool Polymap::getRootVisibilityRegion(const Point2d& source,
                                      VisibilityRegion& visibility_region,
                                      std::string* error) {
  const StopToken no_stop;
  return getRootVisibilityRegion(source, visibility_region, no_stop, error) ==
         OperationStatus::success;
}

OperationStatus Polymap::getRootVisibilityRegion(const Point2d& source,
                                                 VisibilityRegion& visibility_region,
                                                 const StopToken& stop_token,
                                                 std::string* error) {
  visibility_region.clear();
  if (error)
    error->clear();
  if (stop_token.poll())
    return OperationStatus::stopped;

  std::string validation_error;
  const OperationStatus interior_status =
    validateFreeSpaceInterior(source, stop_token, &validation_error);
  if (interior_status != OperationStatus::success) {
    if (interior_status == OperationStatus::stopped)
      return OperationStatus::stopped;
    if (error)
      *error = std::move(validation_error);
    return OperationStatus::failure;
  }

  if (!solution_exist_ || !cdt_ready_) {
    if (error) {
      *error =
        "Root visibility region is unavailable because map geometry "
        "construction failed";
      if (!construction_error_.empty())
        *error += ": " + construction_error_;
    }
    return OperationStatus::failure;
  }
  if (!std::isfinite(source.first) || !std::isfinite(source.second)) {
    if (error)
      *error = "Visibility source is not finite";
    return OperationStatus::failure;
  }

  const exact_geometry::Point exact_source(source.first, source.second);
  // Deliberately pass no obstacle topology.  This both forces closed-cycle
  // validation and bypasses visibility_storage_, whose integer key cannot
  // distinguish two positions in one grid cell.
  if (!calculateVisibilityRegionImpl(
        source, exact_source, std::nullopt, visibility_region, stop_token)) {
    visibility_region.clear();
    if (stop_token.poll()) {
      if (error)
        error->clear();
      return OperationStatus::stopped;
    }
    if (error)
      *error = "Continuous root visibility region calculation failed";
    return OperationStatus::failure;
  }

  if (stop_token.poll()) {
    visibility_region.clear();
    if (error)
      error->clear();
    return OperationStatus::stopped;
  }
  return OperationStatus::success;
}

bool Polymap::validateFreeSpaceInterior(const Point2d& point, std::string* error) const {
  const StopToken no_stop;
  return validateFreeSpaceInterior(point, no_stop, error) == OperationStatus::success;
}

OperationStatus Polymap::validateFreeSpaceInterior(const Point2d& point,
                                                   const StopToken& stop_token,
                                                   std::string* error) const {
  if (error)
    error->clear();
  if (stop_token.poll())
    return OperationStatus::stopped;

  if (validateFreeSpaceInteriorImpl(point, stop_token, error)) {
    if (stop_token.poll()) {
      if (error)
        error->clear();
      return OperationStatus::stopped;
    }
    return OperationStatus::success;
  }
  if (stop_token.poll()) {
    if (error)
      error->clear();
    return OperationStatus::stopped;
  }
  return OperationStatus::failure;
}

bool Polymap::validateFreeSpaceInteriorImpl(const Point2d& point,
                                            const StopToken& stop_token,
                                            std::string* error) const {
  if (error)
    error->clear();
  const auto fail = [&](const std::string& message) {
    if (error)
      *error = message;
    return false;
  };

  if (stop_token.poll())
    return false;
  if (!std::isfinite(point.first) || !std::isfinite(point.second))
    return fail("Free-space endpoint is not finite");
  if (!solution_exist_ || obs_.empty())
    return fail("Free-space boundary is unavailable because the map has no solution");

  const exact_geometry::Point exact_point(point.first, point.second);

  // The flood-fill contour convention keeps the reachable outer contour
  // clockwise (negative signed area) and obstacle holes counter-clockwise.
  // Select the largest-magnitude negative contour, with a largest-area
  // fallback for legacy maps whose winding was not preserved.
  size_t outer_index = 0;
  exact_geometry::FT outer_area(0);
  bool found_clockwise_outer = false;
  for (size_t obstacle_index = 0; obstacle_index < obs_.size(); ++obstacle_index) {
    if (stop_token.poll())
      return false;
    const auto& vertices = obs_[obstacle_index].ordered_vertices_;
    if (vertices.size() < 3)
      return fail("Free-space boundary contains an invalid contour");
    exact_geometry::FT twice_area(0);
    for (size_t index = 0; index < vertices.size(); ++index) {
      if (stop_token.poll())
        return false;
      const auto& from = vertices[index];
      const auto& to = vertices[(index + 1) % vertices.size()];
      twice_area +=
        exact_geometry::FT(from.first) * to.second - exact_geometry::FT(to.first) * from.second;
    }
    if (twice_area == exact_geometry::FT(0))
      return fail("Free-space boundary contains a zero-area contour");
    if ((!found_clockwise_outer && twice_area < exact_geometry::FT(0)) ||
        (found_clockwise_outer && twice_area < outer_area)) {
      outer_index = obstacle_index;
      outer_area = twice_area;
      found_clockwise_outer = true;
    } else if (!found_clockwise_outer && (obstacle_index == 0 || twice_area > outer_area)) {
      outer_index = obstacle_index;
      outer_area = twice_area;
    }
  }

  // Boundary has priority over all point-in-contour classifications so the
  // caller receives an unambiguous boundary diagnostic at hole/outer edges.
  for (size_t obstacle_index = 0; obstacle_index < obs_.size(); ++obstacle_index) {
    if (stop_token.poll())
      return false;
    const auto location =
      locateExactPointInContour(obs_[obstacle_index].ordered_vertices_, exact_point, stop_token);
    if (stop_token.poll())
      return false;
    if (location == ExactContourPointLocation::boundary) {
      return fail("Free-space endpoint lies on an obstacle boundary");
    }
  }

  const auto outer_location =
    locateExactPointInContour(obs_[outer_index].ordered_vertices_, exact_point, stop_token);
  if (stop_token.poll())
    return false;
  if (outer_location != ExactContourPointLocation::inside)
    return fail("Free-space endpoint lies outside the reachable free-space contour");

  for (size_t obstacle_index = 0; obstacle_index < obs_.size(); ++obstacle_index) {
    if (obstacle_index == outer_index)
      continue;
    if (stop_token.poll())
      return false;
    const auto location =
      locateExactPointInContour(obs_[obstacle_index].ordered_vertices_, exact_point, stop_token);
    if (stop_token.poll())
      return false;
    if (location == ExactContourPointLocation::inside)
      return fail("Free-space endpoint lies inside an obstacle contour");
  }
  return true;
}

bool Polymap::getVisibilityRegion(int start_x,
                                  int start_y,
                                  std::vector<std::pair<double, double>>& visibility_region,
                                  std::vector<std::pair<int, int>>& topo_V) {
  return getVisibilityRegion(start_x, start_y, visibility_region, topo_V, nullptr);
}

bool Polymap::getVisibilityRegion(int start_x,
                                  int start_y,
                                  std::vector<std::pair<double, double>>& visibility_region,
                                  std::vector<std::pair<int, int>>& topo_V,
                                  std::string* error) {
  visibility_region.clear();
  topo_V.clear();

  VisibilityRegion rich_region;
  if (!getVisibilityRegion(start_x, start_y, rich_region, error))
    return false;

  visibility_region.reserve(rich_region.size());
  topo_V.reserve(rich_region.size());
  for (const auto& endpoint : rich_region) {
    visibility_region.emplace_back(endpoint.position);
    if (const auto* vertex = std::get_if<ObstacleVertexId>(&endpoint.support))
      topo_V.emplace_back(toLegacyTopology(*vertex));
    else
      topo_V.emplace_back(-1, -1);
  }
  return true;
}

std::pair<int, int> Polymap::locateVertex(int x, int y) const {
  if (xsize_ <= 0 || ysize_ <= 0 || x < 0 || y < 0 || x >= xsize_ || y >= ysize_)
    return {-1, -1};

  const size_t idx = static_cast<size_t>(y) * static_cast<size_t>(xsize_) + static_cast<size_t>(x);
  if (idx >= vertices_location_x_flat_.size() || idx >= vertices_location_y_flat_.size())
    return {-1, -1};

  return {vertices_location_x_flat_[idx], vertices_location_y_flat_[idx]};
}

std::pair<int, int> Polymap::locateVertex(double x, double y) const {
  if (!std::isfinite(x) || !std::isfinite(y))
    return {-1, -1};

  const double rounded_x = std::round(x);
  const double rounded_y = std::round(y);
  constexpr double integer_tolerance = 1e-9;
  if (std::abs(x - rounded_x) > integer_tolerance || std::abs(y - rounded_y) > integer_tolerance ||
      rounded_x < 0.0 || rounded_y < 0.0 || rounded_x >= static_cast<double>(xsize_) ||
      rounded_y >= static_cast<double>(ysize_)) {
    return {-1, -1};
  }

  return locateVertex(static_cast<int>(rounded_x), static_cast<int>(rounded_y));
}

bool Polymap::validateObstacleTopology(std::string& error, bool validate_edge_relations) const {
  const StopToken no_stop;
  return validateObstacleTopologyImpl(error, validate_edge_relations, no_stop);
}

OperationStatus Polymap::validateObstacleTopology(std::string& error,
                                                  bool validate_edge_relations,
                                                  const StopToken& stop_token) const {
  error.clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (validateObstacleTopologyImpl(error, validate_edge_relations, stop_token)) {
    if (stop_token.poll()) {
      error.clear();
      return OperationStatus::stopped;
    }
    return OperationStatus::success;
  }
  if (stop_token.poll()) {
    error.clear();
    return OperationStatus::stopped;
  }
  return OperationStatus::failure;
}

bool Polymap::validateObstacleTopologyImpl(std::string& error,
                                           bool validate_edge_relations,
                                           const StopToken& stop_token) const {
  error.clear();
  const auto fail = [&](std::string message) {
    error = "Unsupported obstacle topology: " + std::move(message);
    return false;
  };

  if (stop_token.poll())
    return false;

  struct ObstacleEdgeRef {
    size_t obstacle;
    size_t edge;
    size_t obstacle_size;
    std::pair<int, int> from_coordinate;
    std::pair<int, int> to_coordinate;
    exact_geometry::Point from;
    exact_geometry::Point to;
  };

  std::unordered_map<uint64_t, ObstacleVertexId> vertex_owners;
  std::vector<ObstacleEdgeRef> edges;
  for (size_t obstacle_index = 0; obstacle_index < obs_.size(); ++obstacle_index) {
    if (stop_token.poll())
      return false;
    const auto& vertices = obs_[obstacle_index].ordered_vertices_;
    if (vertices.size() < 3) {
      return fail("obstacle " + std::to_string(obstacle_index) + " has fewer than three vertices");
    }

    exact_geometry::FT twice_signed_area(0);
    for (size_t vertex_index = 0; vertex_index < vertices.size(); ++vertex_index) {
      if (stop_token.poll())
        return false;
      const auto& vertex = vertices[vertex_index];
      const auto& next = vertices[(vertex_index + 1) % vertices.size()];
      if (vertex.first < 0 || vertex.second < 0 || vertex.first >= xsize_ ||
          vertex.second >= ysize_) {
        return fail("obstacle " + std::to_string(obstacle_index) + " vertex " +
                    std::to_string(vertex_index) + " at (" + std::to_string(vertex.first) + ", " +
                    std::to_string(vertex.second) + ") is outside the vertex registry");
      }

      const uint64_t vertex_key =
        (static_cast<uint64_t>(static_cast<uint32_t>(vertex.first)) << 32) |
        static_cast<uint32_t>(vertex.second);
      const ObstacleVertexId current_owner{static_cast<int>(obstacle_index),
                                           static_cast<int>(vertex_index)};
      const auto inserted = vertex_owners.emplace(vertex_key, current_owner);
      if (!inserted.second) {
        const auto& previous_owner = inserted.first->second;
        if (previous_owner.obstacle != current_owner.obstacle) {
          return fail("obstacles " + std::to_string(previous_owner.obstacle) + " and " +
                      std::to_string(current_owner.obstacle) + " share vertex (" +
                      std::to_string(vertex.first) + ", " + std::to_string(vertex.second) +
                      "); zero-clearance obstacle contacts are not supported");
        }
        return fail("obstacle " + std::to_string(current_owner.obstacle) + " repeats vertex (" +
                    std::to_string(vertex.first) + ", " + std::to_string(vertex.second) +
                    ") at indices " + std::to_string(previous_owner.vertex) + " and " +
                    std::to_string(current_owner.vertex) +
                    "; self-touching contours are not supported");
      }

      if (vertex == next) {
        return fail("obstacle " + std::to_string(obstacle_index) + " edge " +
                    std::to_string(vertex_index) + " has zero length at (" +
                    std::to_string(vertex.first) + ", " + std::to_string(vertex.second) + ")");
      }

      twice_signed_area += exact_geometry::FT(vertex.first) * exact_geometry::FT(next.second) -
                           exact_geometry::FT(next.first) * exact_geometry::FT(vertex.second);
      if (validate_edge_relations) {
        edges.push_back(ObstacleEdgeRef{obstacle_index,
                                        vertex_index,
                                        vertices.size(),
                                        vertex,
                                        next,
                                        exact_geometry::Point(vertex.first, vertex.second),
                                        exact_geometry::Point(next.first, next.second)});
      }
    }

    if (twice_signed_area == exact_geometry::FT(0)) {
      return fail("obstacle " + std::to_string(obstacle_index) + " has zero signed area");
    }
  }

  if (!validate_edge_relations)
    return true;

  for (size_t first_index = 0; first_index < edges.size(); ++first_index) {
    if (stop_token.poll())
      return false;
    const auto& first = edges[first_index];
    for (size_t second_index = first_index + 1; second_index < edges.size(); ++second_index) {
      if (stop_token.poll())
        return false;
      const auto& second = edges[second_index];
      const auto relation = classifyExactSegments(first.from, first.to, second.from, second.to);

      const bool same_obstacle = first.obstacle == second.obstacle;
      const bool first_precedes_second =
        same_obstacle && (first.edge + 1) % first.obstacle_size == second.edge;
      const bool second_precedes_first =
        same_obstacle && (second.edge + 1) % second.obstacle_size == first.edge;
      if (first_precedes_second || second_precedes_first) {
        const exact_geometry::Point& expected_contact =
          first_precedes_second ? first.to : first.from;
        if (relation.relation == ExactSegmentRelation::endpoint_touch && relation.contact &&
            *relation.contact == expected_contact) {
          continue;
        }
      } else if (relation.relation == ExactSegmentRelation::disjoint) {
        continue;
      }

      std::ostringstream message;
      message << "obstacle " << first.obstacle << " edge " << first.edge << " [("
              << first.from_coordinate.first << ", " << first.from_coordinate.second << ") -> ("
              << first.to_coordinate.first << ", " << first.to_coordinate.second
              << ")] and obstacle " << second.obstacle << " edge " << second.edge << " [("
              << second.from_coordinate.first << ", " << second.from_coordinate.second << ") -> ("
              << second.to_coordinate.first << ", " << second.to_coordinate.second << ")] "
              << segmentRelationDescription(relation.relation);
      return fail(message.str());
    }
  }

  return true;
}

bool Polymap::simplificationChordIsTopologicallySafe(size_t obstacle_index,
                                                     size_t current_index) const {
  const StopToken no_stop;
  return simplificationChordIsTopologicallySafeImpl(obstacle_index, current_index, no_stop);
}

OperationStatus Polymap::simplificationChordIsTopologicallySafe(size_t obstacle_index,
                                                                size_t current_index,
                                                                bool& safe,
                                                                const StopToken& stop_token) const {
  safe = false;
  if (stop_token.poll())
    return OperationStatus::stopped;
  safe = simplificationChordIsTopologicallySafeImpl(obstacle_index, current_index, stop_token);
  if (stop_token.poll()) {
    safe = false;
    return OperationStatus::stopped;
  }
  return OperationStatus::success;
}

bool Polymap::simplificationChordIsTopologicallySafeImpl(size_t obstacle_index,
                                                         size_t current_index,
                                                         const StopToken& stop_token) const {
  if (stop_token.poll())
    return false;
  if (obstacle_index >= obs_.size())
    return false;
  const auto& vertices = obs_[obstacle_index].ordered_vertices_;
  if (vertices.size() <= 3 || current_index >= vertices.size())
    return false;

  const size_t previous_index = (current_index + vertices.size() - 1) % vertices.size();
  const size_t next_index = (current_index + 1) % vertices.size();
  const auto& previous_coordinate = vertices[previous_index];
  const auto& next_coordinate = vertices[next_index];
  if (previous_coordinate == next_coordinate)
    return false;

  const exact_geometry::Point previous(previous_coordinate.first, previous_coordinate.second);
  const exact_geometry::Point next(next_coordinate.first, next_coordinate.second);
  const size_t edge_before_previous = (previous_index + vertices.size() - 1) % vertices.size();

  for (size_t other_obstacle = 0; other_obstacle < obs_.size(); ++other_obstacle) {
    if (stop_token.poll())
      return false;
    const auto& other_vertices = obs_[other_obstacle].ordered_vertices_;
    for (size_t edge_index = 0; edge_index < other_vertices.size(); ++edge_index) {
      if (stop_token.poll())
        return false;
      if (other_obstacle == obstacle_index &&
          (edge_index == previous_index || edge_index == current_index)) {
        continue;
      }

      const auto& edge_from_coordinate = other_vertices[edge_index];
      const auto& edge_to_coordinate = other_vertices[(edge_index + 1) % other_vertices.size()];
      if (edge_from_coordinate == edge_to_coordinate)
        return false;
      const exact_geometry::Point edge_from(edge_from_coordinate.first,
                                            edge_from_coordinate.second);
      const exact_geometry::Point edge_to(edge_to_coordinate.first, edge_to_coordinate.second);
      const auto relation = classifyExactSegments(previous, next, edge_from, edge_to);
      if (relation.relation == ExactSegmentRelation::disjoint)
        continue;

      std::optional<exact_geometry::Point> allowed_contact;
      if (other_obstacle == obstacle_index) {
        if (edge_index == edge_before_previous)
          allowed_contact = previous;
        else if (edge_index == next_index)
          allowed_contact = next;
      }
      if (allowed_contact && relation.relation == ExactSegmentRelation::endpoint_touch &&
          relation.contact && *relation.contact == *allowed_contact) {
        continue;
      }
      return false;
    }
  }
  return true;
}

PolymapCreateResult Polymap::create(const GridMap& grid_map,
                                    int start_x,
                                    int start_y,
                                    int goal_x,
                                    int goal_y,
                                    const PlanningLimits& limits) {
  return create(grid_map, start_x, start_y, goal_x, goal_y, StopToken{}, limits);
}

PolymapCreateResult Polymap::create(const GridMap& grid_map,
                                    int start_x,
                                    int start_y,
                                    int goal_x,
                                    int goal_y,
                                    const StopToken& stop_token,
                                    const PlanningLimits& limits) {
  PolymapCreateResult result;
  MapResourceEstimate estimate;
  if (!validateMapResourceBudget(static_cast<size_t>(grid_map.width),
                                 static_cast<size_t>(grid_map.height),
                                 grid_map.data.size(),
                                 limits,
                                 estimate,
                                 result.error)) {
    return result;
  }
  (void)estimate;
  return finishCreation(Polymap(grid_map, start_x, start_y, goal_x, goal_y, stop_token));
}

PolymapCreateResult Polymap::finishCreation(Polymap&& candidate) {
  PolymapCreateResult result;
  if (candidate.construction_stopped_) {
    result.status = PolymapCreateStatus::stopped;
    return result;
  }
  if (candidate.no_path_) {
    result.status = PolymapCreateStatus::no_path;
    result.error = std::move(candidate.construction_error_);
    return result;
  }
  if (!candidate.solution_exist_ || !candidate.cdt_ready_) {
    result.status = PolymapCreateStatus::failure;
    result.error = candidate.construction_error_.empty() ? "Map geometry construction failed"
                                                         : std::move(candidate.construction_error_);
    return result;
  }

  result.status = PolymapCreateStatus::ready;
  result.value.emplace(std::move(candidate));
  return result;
}

Polymap::Polymap(const GridMap& grid_map,
                 int start_x,
                 int start_y,
                 int goal_x,
                 int goal_y,
                 const StopToken& stop_token)
  : xsize_(0), ysize_(0) {
  const auto stop_construction = [&]() { clearStoppedConstructionState(); };

  if (stop_token.poll()) {
    stop_construction();
    return;
  }

  size_t cell_count = 0;
  if (!prepareGridStorage(grid_map, xsize_, ysize_, cell_count, construction_error_)) {
    return;
  }
  data_ = grid_map.data;
  vertices_location_x_flat_.resize(cell_count, -1);
  vertices_location_y_flat_.resize(cell_count, -1);

  if (start_x < 0 || start_y < 0 || goal_x < 0 || goal_y < 0 || start_x >= xsize_ ||
      start_y >= ysize_ || goal_x >= xsize_ || goal_y >= ysize_) {
    construction_error_ = "Start or goal grid cell is outside map bounds";
    return;
  }
  const size_t start_index =
    static_cast<size_t>(start_x) + static_cast<size_t>(start_y) * static_cast<size_t>(xsize_);
  const size_t goal_index =
    static_cast<size_t>(goal_x) + static_cast<size_t>(goal_y) * static_cast<size_t>(xsize_);
  if (data_[start_index] != 0 || data_[goal_index] != 0) {
    construction_error_ = "Start or goal grid cell is occupied";
    return;
  }

  const OperationStatus obstacle_status =
    getPolyObstacles(start_x, start_y, goal_x, goal_y, stop_token);
  if (obstacle_status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  solution_exist_ = obstacle_status == OperationStatus::success;
  if (!solution_exist_) {
    no_path_ = true;
    return;
  }

  OperationStatus status = OperationStatus::success;
  const bool has_outer_contour = hasClockwiseOuterContour(obs_, stop_token);
  if (stop_token.poll()) {
    stop_construction();
    return;
  }

  // Validate the caller's integer endpoints against the raw reachable
  // contours before any vertex simplification can move a boundary edge.  The
  // same check is repeated after simplification below because the simplified
  // contour is the geometry consumed by CDT/visibility.
  std::string endpoint_error;
  if (has_outer_contour) {
    status =
      validateFreeSpaceInterior(Point2d{static_cast<double>(start_x), static_cast<double>(start_y)},
                                stop_token,
                                &endpoint_error);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid start position: " + endpoint_error;
      return;
    }
    status =
      validateFreeSpaceInterior(Point2d{static_cast<double>(goal_x), static_cast<double>(goal_y)},
                                stop_token,
                                &endpoint_error);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid goal position: " + endpoint_error;
      return;
    }
  }

  // Raw contours consist only of unique unit grid edges. An O(V) incidence
  // pass is sufficient here to reject checkerboard/saddle contacts before
  // simplification, without paying an unnecessary O(E^2) exact edge scan.
  status = validateObstacleTopology(construction_error_, false, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;

  status = simplifyPolyObstacles(start_x, start_y, goal_x, goal_y, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure) {
    construction_error_ = "Obstacle simplification failed";
    return;
  }

  status = validateObstacleTopology(construction_error_, true, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;

  if (has_outer_contour) {
    status =
      validateFreeSpaceInterior(Point2d{static_cast<double>(start_x), static_cast<double>(start_y)},
                                stop_token,
                                &endpoint_error);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid start position: " + endpoint_error;
      return;
    }
    status =
      validateFreeSpaceInterior(Point2d{static_cast<double>(goal_x), static_cast<double>(goal_y)},
                                stop_token,
                                &endpoint_error);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid goal position: " + endpoint_error;
      return;
    }
  }

  status = registerVertices(construction_error_, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;

  status = constructCGALRelated(&validateCDTWithCGAL, construction_error_, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;
}

PolymapCreateResult Polymap::create(const GridMap& grid_map,
                                    int start_x,
                                    int start_y,
                                    int goal_x,
                                    int goal_y,
                                    const Point2d& start_position,
                                    const Point2d& goal_position,
                                    const PlanningLimits& limits) {
  return create(
    grid_map, start_x, start_y, goal_x, goal_y, start_position, goal_position, StopToken{}, limits);
}

PolymapCreateResult Polymap::create(const GridMap& grid_map,
                                    int start_x,
                                    int start_y,
                                    int goal_x,
                                    int goal_y,
                                    const Point2d& start_position,
                                    const Point2d& goal_position,
                                    const StopToken& stop_token,
                                    const PlanningLimits& limits) {
  PolymapCreateResult result;
  MapResourceEstimate estimate;
  if (!validateMapResourceBudget(static_cast<size_t>(grid_map.width),
                                 static_cast<size_t>(grid_map.height),
                                 grid_map.data.size(),
                                 limits,
                                 estimate,
                                 result.error)) {
    return result;
  }
  (void)estimate;
  return finishCreation(
    Polymap(grid_map, start_x, start_y, goal_x, goal_y, start_position, goal_position, stop_token));
}

PolymapCreateResult Polymap::create(const GridMap& grid_map,
                                    int start_x,
                                    int start_y,
                                    const Point2d& start_position,
                                    const std::vector<PolymapEndpoint>& goals,
                                    const StopToken& stop_token,
                                    const PlanningLimits& limits) {
  PolymapCreateResult result;
  MapResourceEstimate estimate;
  if (!validateMapResourceBudget(static_cast<size_t>(grid_map.width),
                                 static_cast<size_t>(grid_map.height),
                                 grid_map.data.size(),
                                 limits,
                                 estimate,
                                 result.error)) {
    return result;
  }
  (void)estimate;
  return finishCreation(Polymap(grid_map, start_x, start_y, start_position, goals, stop_token));
}

PolymapCreateResult Polymap::createForReferenceShortening(
  const GridMap& grid_map,
  int start_x,
  int start_y,
  const Point2d& start_position,
  const std::vector<PolymapEndpoint>& goals,
  const StopToken& stop_token,
  const PlanningLimits& limits) {
  PolymapCreateResult result;
  MapResourceEstimate estimate;
  if (!validateMapResourceBudget(static_cast<size_t>(grid_map.width),
                                 static_cast<size_t>(grid_map.height),
                                 grid_map.data.size(),
                                 limits,
                                 estimate,
                                 result.error)) {
    return result;
  }
  const size_t remaining_map_bytes = limits.max_map_bytes - estimate.estimated_bytes;
  const size_t max_raw_contour_vertices =
    remaining_map_bytes / kEstimatedReferenceTopologyBytesPerRawContourVertex;
  return finishCreation(
    Polymap(grid_map,
            start_x,
            start_y,
            start_position,
            goals,
            stop_token,
            false,
            RawContourResourceBudget{max_raw_contour_vertices, limits.max_map_bytes}));
}

Polymap::Polymap(const GridMap& grid_map,
                 int start_x,
                 int start_y,
                 int goal_x,
                 int goal_y,
                 const Point2d& start_position,
                 const Point2d& goal_position,
                 const StopToken& stop_token)
  : Polymap(grid_map,
            start_x,
            start_y,
            start_position,
            std::vector<PolymapEndpoint>{{goal_x, goal_y, goal_position}},
            stop_token) {}

Polymap::Polymap(const GridMap& grid_map,
                 int start_x,
                 int start_y,
                 const Point2d& start_position,
                 const std::vector<PolymapEndpoint>& goals,
                 const StopToken& stop_token)
  : Polymap(grid_map, start_x, start_y, start_position, goals, stop_token, true, std::nullopt) {}

Polymap::Polymap(const GridMap& grid_map,
                 int start_x,
                 int start_y,
                 const Point2d& start_position,
                 const std::vector<PolymapEndpoint>& goals,
                 const StopToken& stop_token,
                 bool simplify_obstacle_contours,
                 std::optional<RawContourResourceBudget> raw_contour_budget)
  : xsize_(0), ysize_(0) {
  const auto stop_construction = [&]() { clearStoppedConstructionState(); };
  const auto reject = [&](const std::string& message) { construction_error_ = message; };

  if (stop_token.poll()) {
    stop_construction();
    return;
  }

  size_t cell_count = 0;
  if (!prepareGridStorage(grid_map, xsize_, ysize_, cell_count, construction_error_)) {
    return;
  }
  data_ = grid_map.data;
  vertices_location_x_flat_.resize(cell_count, -1);
  vertices_location_y_flat_.resize(cell_count, -1);

  if (!std::isfinite(start_position.first) || !std::isfinite(start_position.second)) {
    reject("Start position is not finite");
    return;
  }
  if (goals.empty()) {
    reject("At least one goal is required");
    return;
  }
  if (std::floor(start_position.first) != static_cast<double>(start_x) ||
      std::floor(start_position.second) != static_cast<double>(start_y)) {
    reject("Start continuous position does not belong to the supplied grid cell");
    return;
  }
  if (start_x < 0 || start_y < 0 || start_x >= xsize_ || start_y >= ysize_) {
    reject("Start grid cell is outside map bounds");
    return;
  }
  const size_t start_index =
    static_cast<size_t>(start_x) + static_cast<size_t>(start_y) * static_cast<size_t>(xsize_);
  if (start_index >= data_.size() || data_[start_index] != 0) {
    reject("Start grid cell is occupied");
    return;
  }

  for (const auto& goal : goals) {
    if (!std::isfinite(goal.position.first) || !std::isfinite(goal.position.second)) {
      reject("Goal position is not finite");
      return;
    }
    if (std::floor(goal.position.first) != static_cast<double>(goal.cell_x) ||
        std::floor(goal.position.second) != static_cast<double>(goal.cell_y)) {
      reject("Goal continuous position does not belong to the supplied grid cell");
      return;
    }
    if (goal.cell_x < 0 || goal.cell_y < 0 || goal.cell_x >= xsize_ || goal.cell_y >= ysize_) {
      reject("Goal grid cell is outside map bounds");
      return;
    }
    const size_t goal_index = static_cast<size_t>(goal.cell_x) +
                              static_cast<size_t>(goal.cell_y) * static_cast<size_t>(xsize_);
    if (goal_index >= data_.size() || data_[goal_index] != 0) {
      reject("Goal grid cell is occupied");
      return;
    }
  }

  const std::optional<size_t> max_raw_contour_vertices =
    raw_contour_budget ? std::optional<size_t>(raw_contour_budget->max_vertices) : std::nullopt;
  const OperationStatus obstacle_status =
    getPolyObstacles(start_x, start_y, goals, stop_token, max_raw_contour_vertices);
  if (obstacle_status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  solution_exist_ = obstacle_status == OperationStatus::success;
  if (!solution_exist_) {
    if (!construction_error_.empty()) {
      if (raw_contour_budget) {
        construction_error_ +=
          " derived from max_map_bytes=" +
          std::to_string(raw_contour_budget->max_map_bytes);
      }
      return;
    }
    no_path_ = true;
    construction_error_ = "Start and every goal must be in the same reachable free-space component";
    return;
  }

  if (!simplify_obstacle_contours && raw_contour_budget) {
    const size_t max_raw_contour_vertices = raw_contour_budget->max_vertices;
    size_t raw_contour_vertices = 0;
    for (const auto& obstacle : obs_) {
      if (stop_token.poll()) {
        stop_construction();
        return;
      }
      if (raw_contour_vertices > max_raw_contour_vertices ||
          obstacle.ordered_vertices_.size() >
            max_raw_contour_vertices - raw_contour_vertices) {
        construction_error_ =
          "Reference-shortening raw contour exceeds the vertex budget " +
          std::to_string(max_raw_contour_vertices) + " derived from max_map_bytes=" +
          std::to_string(raw_contour_budget->max_map_bytes);
        return;
      }
      raw_contour_vertices += obstacle.ordered_vertices_.size();
    }
  }

  OperationStatus status = OperationStatus::success;
  const bool has_outer_contour = hasClockwiseOuterContour(obs_, stop_token);
  if (stop_token.poll()) {
    stop_construction();
    return;
  }

  // Reject continuous endpoints against the unsimplified contours first.  In
  // particular, a point on a raw contour must never become accepted merely
  // because a later chord removes the incident contour vertices.
  std::string endpoint_error;
  if (has_outer_contour) {
    status = validateFreeSpaceInterior(start_position, stop_token, &endpoint_error);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid start position: " + endpoint_error;
      return;
    }
    for (const auto& goal : goals) {
      status = validateFreeSpaceInterior(goal.position, stop_token, &endpoint_error);
      if (status == OperationStatus::stopped) {
        stop_construction();
        return;
      }
      if (status == OperationStatus::failure) {
        construction_error_ = "Invalid goal position: " + endpoint_error;
        return;
      }
    }
  }

  status = validateObstacleTopology(construction_error_, false, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;

  if (simplify_obstacle_contours) {
    std::vector<Point2d> protected_points;
    protected_points.reserve(goals.size() + 1);
    protected_points.emplace_back(start_position);
    for (const auto& goal : goals) protected_points.emplace_back(goal.position);
    status = simplifyPolyObstaclesImpl(protected_points, stop_token)
               ? OperationStatus::success
               : (stop_token.poll() ? OperationStatus::stopped : OperationStatus::failure);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Obstacle simplification failed";
      return;
    }

    status = validateObstacleTopology(construction_error_, true, stop_token);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure)
      return;

    if (has_outer_contour) {
      status = validateFreeSpaceInterior(start_position, stop_token, &endpoint_error);
      if (status == OperationStatus::stopped) {
        stop_construction();
        return;
      }
      if (status == OperationStatus::failure) {
        construction_error_ = "Invalid start position: " + endpoint_error;
        return;
      }
      for (const auto& goal : goals) {
        status = validateFreeSpaceInterior(goal.position, stop_token, &endpoint_error);
        if (status == OperationStatus::stopped) {
          stop_construction();
          return;
        }
        if (status == OperationStatus::failure) {
          construction_error_ = "Invalid goal position: " + endpoint_error;
          return;
        }
      }
    }
  }

  status = registerVertices(construction_error_, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;

  status = constructCGALRelated(&validateCDTWithCGAL, construction_error_, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;
}

bool Polymap::registerVertices(std::string& error) {
  const StopToken no_stop;
  return registerVerticesImpl(error, no_stop);
}

OperationStatus Polymap::registerVertices(std::string& error, const StopToken& stop_token) {
  error.clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (registerVerticesImpl(error, stop_token))
    return OperationStatus::success;
  if (stop_token.poll()) {
    error.clear();
    return OperationStatus::stopped;
  }
  return OperationStatus::failure;
}

bool Polymap::registerVerticesImpl(std::string& error, const StopToken& stop_token) {
  error.clear();
  if (stop_token.poll())
    return false;
  const size_t registry_size = static_cast<size_t>(xsize_) * static_cast<size_t>(ysize_);
  std::vector<int> candidate_obstacles(registry_size, -1);
  std::vector<int> candidate_vertices(registry_size, -1);

  for (size_t i = 0; i < obs_.size(); ++i) {
    if (stop_token.poll())
      return false;
    for (size_t j = 0; j < obs_[i].ordered_vertices_.size(); ++j) {
      if (stop_token.poll())
        return false;
      const int x = obs_[i].ordered_vertices_[j].first;
      const int y = obs_[i].ordered_vertices_[j].second;
      if (x < 0 || y < 0 || x >= xsize_ || y >= ysize_) {
        error = "Unsupported obstacle topology: obstacle " + std::to_string(i) + " vertex " +
                std::to_string(j) + " at (" + std::to_string(x) + ", " + std::to_string(y) +
                ") is outside the vertex registry";
        return false;
      }

      const size_t index =
        static_cast<size_t>(y) * static_cast<size_t>(xsize_) + static_cast<size_t>(x);
      if (candidate_obstacles[index] >= 0) {
        error = "Unsupported obstacle topology: obstacle " + std::to_string(i) + " vertex " +
                std::to_string(j) + " would overwrite obstacle " +
                std::to_string(candidate_obstacles[index]) + " vertex " +
                std::to_string(candidate_vertices[index]) + " at shared vertex (" +
                std::to_string(x) + ", " + std::to_string(y) + ")";
        return false;
      }
      candidate_obstacles[index] = static_cast<int>(i);
      candidate_vertices[index] = static_cast<int>(j);
    }
  }

  if (stop_token.poll())
    return false;
  vertices_location_x_flat_.swap(candidate_obstacles);
  vertices_location_y_flat_.swap(candidate_vertices);
  return true;
}

void Polymap::simplifyPolyObstacles(int start_x, int start_y, int goal_x, int goal_y) {
  const StopToken no_stop;
  (void)simplifyPolyObstaclesImpl(
    Point2d{static_cast<double>(start_x), static_cast<double>(start_y)},
    Point2d{static_cast<double>(goal_x), static_cast<double>(goal_y)},
    no_stop);
}

OperationStatus Polymap::simplifyPolyObstacles(
  int start_x, int start_y, int goal_x, int goal_y, const StopToken& stop_token) {
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (simplifyPolyObstaclesImpl(Point2d{static_cast<double>(start_x), static_cast<double>(start_y)},
                                Point2d{static_cast<double>(goal_x), static_cast<double>(goal_y)},
                                stop_token)) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    return OperationStatus::success;
  }
  return stop_token.poll() ? OperationStatus::stopped : OperationStatus::failure;
}

OperationStatus Polymap::simplifyPolyObstacles(const Point2d& start,
                                               const Point2d& goal,
                                               const StopToken& stop_token) {
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (simplifyPolyObstaclesImpl(start, goal, stop_token)) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    return OperationStatus::success;
  }
  return stop_token.poll() ? OperationStatus::stopped : OperationStatus::failure;
}

void Polymap::simplifyPolyObstacles(const Point2d& start, const Point2d& goal) {
  const StopToken no_stop;
  (void)simplifyPolyObstaclesImpl(start, goal, no_stop);
}

bool Polymap::simplifyPolyObstaclesImpl(const Point2d& start,
                                        const Point2d& goal,
                                        const StopToken& stop_token) {
  return simplifyPolyObstaclesImpl(std::vector<Point2d>{start, goal}, stop_token);
}

bool Polymap::simplifyPolyObstaclesImpl(const std::vector<Point2d>& protected_points,
                                        const StopToken& stop_token) {
  for (auto iter = obs_.begin(); iter != obs_.end(); ++iter) {
    if (stop_token.poll())
      return false;
    int prev, curr, next;
    bool stable = false;
    curr = 0;
    bool simplifable;
    int x1, y1, x2, y2, x3, y3;

    while (1) {
      if (stop_token.poll())
        return false;
      if (iter->ordered_vertices_.size() <= 3)
        break;
      prev = (curr - 1 + static_cast<int>(iter->ordered_vertices_.size())) %
             static_cast<int>(iter->ordered_vertices_.size());
      next = (curr + 1) % iter->ordered_vertices_.size();

      x1 = iter->ordered_vertices_[prev].first;
      y1 = iter->ordered_vertices_[prev].second;
      x2 = iter->ordered_vertices_[curr].first;
      y2 = iter->ordered_vertices_[curr].second;
      x3 = iter->ordered_vertices_[next].first;
      y3 = iter->ordered_vertices_[next].second;

      const exact_geometry::Point previous_point(x1, y1);
      const exact_geometry::Point current_point(x2, y2);
      const exact_geometry::Point next_point(x3, y3);
      const CGAL::Orientation turn = CGAL::orientation(previous_point, current_point, next_point);

      if (turn == CGAL::COLLINEAR) {
        // Only a true middle point is redundant. A collinear backtracking
        // spike changes the contour and must not be erased.
        simplifable =
          exact_geometry::isRemovableCollinearMiddle(previous_point, current_point, next_point);
      } else {
        simplifable = true;
        if (turn == CGAL::RIGHT_TURN) {
          for (const auto& point : protected_points) {
            if (isInTri(x1, y1, x2, y2, x3, y3, point.first, point.second)) {
              simplifable = false;
              break;
            }
          }

          double testx, testy;
          if (simplifable) {
            for (auto iter2 = obs_.begin(); iter2 != obs_.end(); ++iter2) {
              if (stop_token.poll())
                return false;
              for (auto iter3 = iter2->ordered_vertices_.begin();
                   iter3 != iter2->ordered_vertices_.end();
                   ++iter3) {
                if (stop_token.poll())
                  return false;
                testx = iter3->first;
                testy = iter3->second;
                if (isInTri(x1, y1, x2, y2, x3, y3, testx, testy)) {
                  if (iter2 - obs_.begin() != iter - obs_.begin())
                    simplifable = false;
                  else {
                    if (iter3 - iter2->ordered_vertices_.begin() != prev &&
                        iter3 - iter2->ordered_vertices_.begin() != curr &&
                        iter3 - iter2->ordered_vertices_.begin() != next)
                      simplifable = false;
                  }
                }
              }
            }
          }
        } else {
          simplifable = false;
        }
      }

      if (simplifable) {
        const size_t obstacle_index = static_cast<size_t>(iter - obs_.begin());
        simplifable = simplificationChordIsTopologicallySafeImpl(
          obstacle_index, static_cast<size_t>(curr), stop_token);
        if (stop_token.poll())
          return false;
      }

      if (simplifable) {
        iter->ordered_vertices_.erase(iter->ordered_vertices_.begin() + curr);
        if (curr >= static_cast<int>(iter->ordered_vertices_.size()))
          curr = static_cast<int>(iter->ordered_vertices_.size()) - 1;
        stable = false;
        continue;
      } else {
        if (curr == 0) {
          if (!stable)
            stable = true;
          else
            break;
        }
      }

      curr++;
      if (curr >= static_cast<int>(iter->ordered_vertices_.size()))
        curr = 0;
    }
  }
  return true;
}

OperationStatus Polymap::buildTriangleEnvironment(std::string& error,
                                                  const StopToken& stop_token) {
  using EdgeKey = std::pair<Point2d, Point2d>;
  struct EdgeUse {
    int face = -1;
    int local_edge = -1;
  };
  struct EdgeAccumulator {
    std::vector<EdgeUse> uses;
    bool constrained = false;
  };

  error.clear();
  triangle_faces_.clear();
  triangle_edges_.clear();
  const auto stopped = [&]() {
    error.clear();
    triangle_faces_.clear();
    triangle_edges_.clear();
    return OperationStatus::stopped;
  };
  const auto fail = [&](std::string message) {
    error = std::move(message);
    triangle_faces_.clear();
    triangle_edges_.clear();
    return OperationStatus::failure;
  };
  if (stop_token.poll())
    return stopped();
  if (!cdt_ready_ || facets_.empty() || obs_.empty()) {
    return fail("Triangle environment requires a ready CDT and a reachable outer contour");
  }
  if (facets_.size() > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
    return fail("Triangle count exceeds the stable TriangleId range");
  }

  const auto canonical_edge = [](Point2d first, Point2d second) {
    if (second < first)
      std::swap(first, second);
    return EdgeKey{first, second};
  };

  std::set<EdgeKey> constrained_edges;
  for (const auto& obstacle : obs_) {
    if (stop_token.poll())
      return stopped();
    const auto& vertices = obstacle.ordered_vertices_;
    for (size_t index = 0; index < vertices.size(); ++index) {
      if (stop_token.poll())
        return stopped();
      const auto& first = vertices[index];
      const auto& second = vertices[(index + 1) % vertices.size()];
      constrained_edges.insert(
        canonical_edge({static_cast<double>(first.first), static_cast<double>(first.second)},
                       {static_cast<double>(second.first), static_cast<double>(second.second)}));
    }
  }

  size_t outer_index = 0;
  long double outer_area = 0.0L;
  bool found_clockwise_outer = false;
  for (size_t obstacle_index = 0; obstacle_index < obs_.size(); ++obstacle_index) {
    if (stop_token.poll())
      return stopped();
    const auto& vertices = obs_[obstacle_index].ordered_vertices_;
    long double twice_area = 0.0L;
    for (size_t index = 0; index < vertices.size(); ++index) {
      if (stop_token.poll())
        return stopped();
      const auto& from = vertices[index];
      const auto& to = vertices[(index + 1) % vertices.size()];
      twice_area += static_cast<long double>(from.first) * to.second -
                    static_cast<long double>(to.first) * from.second;
    }
    if ((!found_clockwise_outer && twice_area < 0.0L) ||
        (found_clockwise_outer && twice_area < outer_area)) {
      outer_index = obstacle_index;
      outer_area = twice_area;
      found_clockwise_outer = true;
    } else if (!found_clockwise_outer && (obstacle_index == 0 || twice_area > outer_area)) {
      outer_index = obstacle_index;
      outer_area = twice_area;
    }
  }

  triangle_faces_.resize(facets_.size());
  std::map<EdgeKey, EdgeAccumulator> edges;
  for (size_t face_index = 0; face_index < facets_.size(); ++face_index) {
    if (stop_token.poll())
      return stopped();
    const auto& facet = facets_[face_index];
    if (facet.size() != 3) {
      return fail("CDT facet does not contain exactly three vertices");
    }
    auto& face = triangle_faces_[face_index];
    long double centroid_x = 0.0L;
    long double centroid_y = 0.0L;
    for (size_t vertex_index = 0; vertex_index < 3; ++vertex_index) {
      face.vertices[vertex_index] = {static_cast<double>(facet[vertex_index].first),
                                     static_cast<double>(facet[vertex_index].second)};
      centroid_x += facet[vertex_index].first;
      centroid_y += facet[vertex_index].second;
    }
    centroid_x /= 3.0L;
    centroid_y /= 3.0L;
    bool point_is_inside = false;
    if (pointInPolygon(obs_[outer_index].ordered_vertices_,
                       centroid_x,
                       centroid_y,
                       stop_token,
                       point_is_inside) == OperationStatus::stopped) {
      return stopped();
    }
    face.is_free = point_is_inside;
    for (size_t obstacle_index = 0; obstacle_index < obs_.size() && face.is_free;
         ++obstacle_index) {
      if (stop_token.poll())
        return stopped();
      if (obstacle_index == outer_index)
        continue;
      if (pointInPolygon(obs_[obstacle_index].ordered_vertices_,
                         centroid_x,
                         centroid_y,
                         stop_token,
                         point_is_inside) == OperationStatus::stopped) {
        return stopped();
      }
      if (point_is_inside)
        face.is_free = false;
    }

    for (size_t edge_index = 0; edge_index < 3; ++edge_index) {
      if (stop_token.poll())
        return stopped();
      const EdgeKey key =
        canonical_edge(face.vertices[edge_index], face.vertices[(edge_index + 1) % 3]);
      auto& accumulator = edges[key];
      accumulator.uses.push_back({static_cast<int>(face_index), static_cast<int>(edge_index)});
      accumulator.constrained = constrained_edges.find(key) != constrained_edges.end();
    }
  }

  for (auto& [key, accumulator] : edges) {
    if (stop_token.poll())
      return stopped();
    if (accumulator.uses.empty() || accumulator.uses.size() > 2) {
      return fail("CDT edge has an invalid number of incident finite faces");
    }
    TriangleMeshEdge edge;
    edge.a = key.first;
    edge.b = key.second;
    edge.constrained = accumulator.constrained;
    for (size_t use_index = 0; use_index < accumulator.uses.size(); ++use_index) {
      const auto& use = accumulator.uses[use_index];
      edge.faces[use_index] = use.face;
      triangle_faces_[static_cast<size_t>(use.face)]
        .constrained[static_cast<size_t>(use.local_edge)] = edge.constrained;
    }
    if (accumulator.uses.size() == 2 && !edge.constrained) {
      const auto& first = accumulator.uses[0];
      const auto& second = accumulator.uses[1];
      if (triangle_faces_[static_cast<size_t>(first.face)].is_free &&
          triangle_faces_[static_cast<size_t>(second.face)].is_free) {
        triangle_faces_[static_cast<size_t>(first.face)]
          .neighbors[static_cast<size_t>(first.local_edge)] = second.face;
        triangle_faces_[static_cast<size_t>(second.face)]
          .neighbors[static_cast<size_t>(second.local_edge)] = first.face;
      }
    }
    triangle_edges_.push_back(edge);
  }

  size_t free_count = 0;
  for (auto fit = cdt_.finite_faces_begin(); fit != cdt_.finite_faces_end(); ++fit) {
    if (stop_token.poll())
      return stopped();
    const int id = fit->info().stable_id;
    if (id < 0 || id >= static_cast<int>(triangle_faces_.size())) {
      return fail("CDT face lost its stable ID while building the triangle environment");
    }
    fit->info().is_free = triangle_faces_[static_cast<size_t>(id)].is_free;
    if (fit->info().is_free)
      ++free_count;
  }
  if (free_count == 0) {
    return fail("Triangle environment contains no free faces");
  }
  if (stop_token.poll())
    return stopped();
  return OperationStatus::success;
}

size_t Polymap::freeTriangleCount() const noexcept {
  return static_cast<size_t>(std::count_if(
    triangle_faces_.begin(), triangle_faces_.end(), [](const auto& face) { return face.is_free; }));
}

OperationStatus Polymap::traceFreeSpacePath(const std::vector<Point2d>& path,
                                            TriangleCorridor& corridor,
                                            const StopToken& stop_token,
                                            std::string* error) const {
  corridor = {};
  if (error)
    error->clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  return traceFreeSpacePathImpl(path, corridor, stop_token, error);
}

OperationStatus Polymap::traceFreeSpacePathImpl(const std::vector<Point2d>& input_path,
                                                TriangleCorridor& corridor,
                                                const StopToken& stop_token,
                                                std::string* error) const {
  const auto stopped = [&]() {
    corridor = {};
    if (error)
      error->clear();
    return OperationStatus::stopped;
  };
  const auto fail = [&](const std::string& message) {
    corridor = {};
    if (error)
      *error = message;
    return OperationStatus::failure;
  };
  if (!cdt_ready_ || triangle_faces_.empty())
    return fail("Free-triangle environment is not ready");
  if (input_path.empty())
    return fail("A reference polyline requires at least one point");

  std::vector<Point2d> path;
  path.reserve(input_path.size());
  for (const auto& point : input_path) {
    if (!std::isfinite(point.first) || !std::isfinite(point.second))
      return fail("Reference polyline contains a non-finite point");
    if (path.empty() || !samePoint(path.back(), point))
      path.push_back(point);
  }
  const auto locate_candidates = [&](const Point2d& point) {
    std::vector<int> candidates;
    using CDT = constrained_delaunay_triangulation::CDT;
    CDT::Locate_type locate_type;
    int local_index = -1;
    const auto face =
      cdt_.locate(constrained_delaunay_triangulation::Point(point.first, point.second),
                  locate_type,
                  local_index);
    const auto add_face = [&](const CDT::Face_handle& candidate) {
      if (candidate == CDT::Face_handle() || cdt_.is_infinite(candidate))
        return;
      const int id = candidate->info().stable_id;
      if (id >= 0 && id < static_cast<int>(triangle_faces_.size()) &&
          triangle_faces_[static_cast<size_t>(id)].is_free) {
        candidates.push_back(id);
      }
    };
    if (locate_type == CDT::FACE) {
      add_face(face);
    } else if (locate_type == CDT::EDGE) {
      add_face(face);
      add_face(face->neighbor(local_index));
    } else if (locate_type == CDT::VERTEX) {
      const auto vertex = face->vertex(local_index);
      auto incident = cdt_.incident_faces(vertex, face);
      const auto done = incident;
      if (incident != 0) {
        do {
          add_face(incident);
          ++incident;
        } while (incident != done);
      }
    }
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    return candidates;
  };

  // A constant path is the identity transition. Keep one containing triangle
  // as its zero-portal sleeve so callers can distinguish it from an invalid
  // or out-of-free-space reference.
  if (path.size() == 1) {
    const auto candidates = locate_candidates(path.front());
    if (candidates.empty())
      return fail("The constant reference point is outside the free-triangle environment");
    corridor.triangle_occurrences.push_back(static_cast<uint32_t>(candidates.front()));
    return OperationStatus::success;
  }

  const auto connect_at_point = [&](int from, int to, const Point2d& point) {
    std::vector<int> empty;
    if (from == to)
      return std::vector<int>{from};
    std::vector<int> parent(triangle_faces_.size(), -2);
    std::queue<int> pending;
    parent[static_cast<size_t>(from)] = -1;
    pending.push(from);
    while (!pending.empty()) {
      if (stop_token.poll())
        return empty;
      const int current = pending.front();
      pending.pop();
      std::vector<int> neighbors;
      for (const int neighbor : triangle_faces_[static_cast<size_t>(current)].neighbors) {
        if (neighbor >= 0)
          neighbors.push_back(neighbor);
      }
      std::sort(neighbors.begin(), neighbors.end());
      neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
      for (const int neighbor : neighbors) {
        if (parent[static_cast<size_t>(neighbor)] != -2 ||
            !pointInTriangleClosed(triangle_faces_[static_cast<size_t>(neighbor)].vertices,
                                   point)) {
          continue;
        }
        parent[static_cast<size_t>(neighbor)] = current;
        if (neighbor == to) {
          std::vector<int> connection;
          for (int trace = to; trace >= 0; trace = parent[static_cast<size_t>(trace)])
            connection.push_back(trace);
          std::reverse(connection.begin(), connection.end());
          return connection;
        }
        pending.push(neighbor);
      }
    }
    return empty;
  };

  std::string segment_failure;
  const auto segment_faces = [&](const Point2d& start, const Point2d& goal) {
    segment_failure.clear();
    std::vector<int> faces;
    std::vector<long double> events{0.0L, 1.0L};
    events.reserve(triangle_edges_.size() + 2);
    for (const auto& edge : triangle_edges_) {
      if (stop_token.poll())
        return faces;
      appendSegmentEdgeEvents(start, goal, edge.a, edge.b, events);
    }
    std::sort(events.begin(), events.end());
    constexpr long double event_tolerance = 1.0e-13L;
    events.erase(std::unique(events.begin(),
                             events.end(),
                             [](long double first, long double second) {
                               return std::abs(first - second) <= event_tolerance;
                             }),
                 events.end());

    for (size_t interval = 0; interval + 1 < events.size(); ++interval) {
      if (stop_token.poll())
        return std::vector<int>{};
      if (events[interval + 1] - events[interval] <= event_tolerance)
        continue;
      const long double parameter = (events[interval] + events[interval + 1]) / 2.0L;
      const Point2d midpoint{
        static_cast<double>(static_cast<long double>(start.first) +
                            parameter * (static_cast<long double>(goal.first) - start.first)),
        static_cast<double>(static_cast<long double>(start.second) +
                            parameter * (static_cast<long double>(goal.second) - start.second))};
      const auto candidates = locate_candidates(midpoint);
      if (candidates.empty()) {
        std::ostringstream stream;
        stream << "no free face contains interval midpoint (" << midpoint.first << ", "
               << midpoint.second << ")";
        segment_failure = stream.str();
        return std::vector<int>{};
      }

      int selected = candidates.front();
      std::vector<int> selected_connection;
      if (!faces.empty()) {
        const Point2d event_point{
          static_cast<double>(static_cast<long double>(start.first) +
                              events[interval] *
                                (static_cast<long double>(goal.first) - start.first)),
          static_cast<double>(static_cast<long double>(start.second) +
                              events[interval] *
                                (static_cast<long double>(goal.second) - start.second))};
        size_t best_size = std::numeric_limits<size_t>::max();
        for (const int candidate : candidates) {
          auto connection = connect_at_point(faces.back(), candidate, event_point);
          if (!connection.empty() && (connection.size() < best_size ||
                                      (connection.size() == best_size && candidate < selected))) {
            selected = candidate;
            best_size = connection.size();
            selected_connection = std::move(connection);
          }
        }
        if (selected_connection.empty()) {
          std::ostringstream stream;
          stream << "no incident free-face connection from triangle " << faces.back()
                 << " at parameter " << static_cast<double>(events[interval]);
          segment_failure = stream.str();
          return std::vector<int>{};
        }
        faces.insert(faces.end(), selected_connection.begin() + 1, selected_connection.end());
      }
      if (faces.empty() || faces.back() != selected)
        faces.push_back(selected);
    }
    return faces;
  };

  std::vector<int> raw_faces;
  for (size_t segment_index = 1; segment_index < path.size(); ++segment_index) {
    if (stop_token.poll())
      return stopped();
    auto current_faces = segment_faces(path[segment_index - 1], path[segment_index]);
    // segment_faces() and its connect_at_point() helper use an empty vector
    // for both geometric failure and cooperative stop.  Check the latched
    // token before interpreting that sentinel as a missing free-space sleeve.
    if (stop_token.poll())
      return stopped();
    if (current_faces.empty()) {
      return fail("Reference segment " + std::to_string(segment_index - 1) + "->" +
                  std::to_string(segment_index) + " leaves the free-triangle environment" +
                  (segment_failure.empty() ? std::string{} : ": " + segment_failure));
    }
    if (!raw_faces.empty()) {
      auto connection =
        connect_at_point(raw_faces.back(), current_faces.front(), path[segment_index - 1]);
      if (stop_token.poll())
        return stopped();
      if (connection.empty())
        return fail("Consecutive reference segments do not share a free triangle fan");
      raw_faces.insert(raw_faces.end(), connection.begin() + 1, connection.end());
    }
    if (raw_faces.empty()) {
      raw_faces.insert(raw_faces.end(), current_faces.begin(), current_faces.end());
    } else {
      const size_t offset = raw_faces.back() == current_faces.front() ? 1 : 0;
      raw_faces.insert(raw_faces.end(), current_faces.begin() + offset, current_faces.end());
    }
  }

  std::vector<int> reduced_faces;
  reduced_faces.reserve(raw_faces.size());
  for (const int face : raw_faces) {
    if (!reduced_faces.empty() && reduced_faces.back() == face)
      continue;
    if (reduced_faces.size() >= 2 && reduced_faces[reduced_faces.size() - 2] == face) {
      reduced_faces.pop_back();
      continue;
    }
    reduced_faces.push_back(face);
  }
  if (reduced_faces.empty())
    return fail("Reference path produced an empty triangle corridor");

  corridor.triangle_occurrences.reserve(reduced_faces.size());
  for (const int face : reduced_faces)
    corridor.triangle_occurrences.push_back(static_cast<uint32_t>(face));
  corridor.portals.reserve(reduced_faces.size() - 1);
  for (size_t index = 1; index < reduced_faces.size(); ++index) {
    const int from = reduced_faces[index - 1];
    const int to = reduced_faces[index];
    const auto& from_face = triangle_faces_[static_cast<size_t>(from)];
    int shared_edge = -1;
    for (size_t edge_index = 0; edge_index < 3; ++edge_index) {
      if (from_face.neighbors[edge_index] == to) {
        shared_edge = static_cast<int>(edge_index);
        break;
      }
    }
    if (shared_edge < 0)
      return fail("Reduced triangle corridor contains non-adjacent faces");
    const Point2d first = from_face.vertices[static_cast<size_t>(shared_edge)];
    const Point2d second = from_face.vertices[(static_cast<size_t>(shared_edge) + 1) % 3];
    DirectedPortal portal;
    portal.from_triangle = static_cast<uint32_t>(from);
    portal.to_triangle = static_cast<uint32_t>(to);
    // The common edge is stored in the boundary order of from_face. Unlike
    // a centroid-to-centroid test, this remains well defined for obtuse
    // adjacent triangles where both portal endpoints can lie on the same
    // side of the centroid line. For a counter-clockwise source face the
    // source interior is left of first->second, so crossing to the neighbor
    // makes second the portal's geometric left endpoint and first its right.
    // Reversing the traversal therefore swaps the endpoints exactly.
    const int face_orientation =
      exactOrientation(from_face.vertices[0], from_face.vertices[1], from_face.vertices[2]);
    if (face_orientation > 0) {
      portal.left = second;
      portal.right = first;
    } else {
      portal.left = first;
      portal.right = second;
    }
    corridor.portals.push_back(portal);
  }
  return OperationStatus::success;
}

HomotopyShorteningResult Polymap::shortenPathWithinHomotopy(const std::vector<Point2d>& reference,
                                                            const StopToken& stop_token) const {
  HomotopyShorteningResult result;
  const auto clear_uncertified_output = [&result]() {
    result.path.clear();
    result.path_cost = 0.0;
    result.collision_free = false;
    result.homotopy_preserved = false;
    result.locally_shortest = false;
  };
  if (stop_token.poll()) {
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled before tracing";
    return result;
  }
  std::string trace_error;
  const OperationStatus trace_status =
    traceFreeSpacePath(reference, result.corridor, stop_token, &trace_error);
  if (trace_status == OperationStatus::stopped) {
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled while tracing the reference";
    return result;
  }
  if (trace_status != OperationStatus::success) {
    result.status = reference.empty() ? HomotopyShorteningStatus::invalid_reference
                                      : HomotopyShorteningStatus::no_corridor;
    result.message = trace_error.empty() ? "Could not trace the reference path" : trace_error;
    return result;
  }
  std::string witness_error;
  if (!validateReducedDirectedPortalWitness(result.corridor, &witness_error)) {
    result.status = HomotopyShorteningStatus::failure;
    result.message = "Reference trace produced an invalid lifted-sleeve witness: " + witness_error;
    return result;
  }
  if (stop_token.poll()) {
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled before funnel execution";
    return result;
  }

  result.path = runFunnel(reference.front(), reference.back(), result.corridor.portals);
  if (result.path.empty() || !samePoint(result.path.front(), reference.front()) ||
      !samePoint(result.path.back(), reference.back())) {
    clear_uncertified_output();
    result.status = HomotopyShorteningStatus::failure;
    result.message = "Funnel result did not preserve the exact reference endpoints";
    return result;
  }
  result.path_cost = polylineLength(result.path);
  const double reference_cost = polylineLength(reference);
  const double cost_tolerance = 1.0e-10 * std::max({1.0, reference_cost, result.path_cost});
  if (result.path_cost > reference_cost + cost_tolerance) {
    std::ostringstream diagnostic;
    diagnostic << "Funnel result cost " << result.path_cost << " exceeds reference cost "
               << reference_cost << "; output=";
    for (const auto& point : result.path)
      diagnostic << "(" << point.first << "," << point.second << ")";
    diagnostic << "; portals=";
    for (const auto& portal : result.corridor.portals) {
      diagnostic << "[" << portal.from_triangle << "->" << portal.to_triangle << " L("
                 << portal.left.first << "," << portal.left.second << ") R(" << portal.right.first
                 << "," << portal.right.second << ")]";
    }
    clear_uncertified_output();
    result.status = HomotopyShorteningStatus::failure;
    result.message = diagnostic.str();
    return result;
  }

  std::string output_error;
  const OperationStatus output_status =
    traceFreeSpacePath(result.path, result.output_corridor, stop_token, &output_error);
  if (output_status == OperationStatus::stopped) {
    clear_uncertified_output();
    result.output_corridor = {};
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled during output validation";
    return result;
  }
  if (output_status != OperationStatus::success) {
    clear_uncertified_output();
    result.output_corridor = {};
    result.status = HomotopyShorteningStatus::failure;
    result.message = "Funnel result failed free-space validation: " + output_error;
    return result;
  }
  if (stop_token.poll()) {
    clear_uncertified_output();
    result.output_corridor = {};
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled before portal-witness certification";
    return result;
  }
  // The contour CDT adds no free-space Steiner vertices: every free face is
  // convex and every unconstrained shared edge is a contractible portal.
  // Contracting path fragments inside successive faces therefore maps a path
  // to its reduced ordered dual-edge walk.  Equality of the complete directed
  // portal-occurrence walk (not merely a set of face IDs) is a sufficient
  // relative-endpoint homotopy certificate, including repeated winding
  // cycles.  Re-tracing the emitted geometry makes this proof independent of
  // the funnel implementation itself.
  if (!sameReducedDirectedPortalWitness(result.corridor, result.output_corridor, &witness_error)) {
    clear_uncertified_output();
    result.status = HomotopyShorteningStatus::failure;
    result.message = "Funnel result changed the reduced lifted-sleeve witness: " + witness_error;
    return result;
  }
  if (stop_token.poll()) {
    clear_uncertified_output();
    result.output_corridor = {};
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled after portal-witness certification";
    return result;
  }

  result.status = HomotopyShorteningStatus::success;
  result.collision_free = true;
  result.homotopy_preserved = true;
  result.locally_shortest = true;
  result.message =
    "Funnel completed and independently retraced the same reduced directed-portal witness";
  return result;
}

bool Polymap::isFacetInsideObstacle(int facet_idx) const {
  if (facet_idx < 0 || facet_idx >= static_cast<int>(facets_.size()))
    return false;
  const auto& f = facets_[facet_idx];
  double cx = (f[0].first + f[1].first + f[2].first) / 3.0;
  double cy = (f[0].second + f[1].second + f[2].second) / 3.0;
  for (size_t oi = 1; oi < obs_.size(); ++oi) {
    const auto& v = obs_[oi].ordered_vertices_;
    int n = static_cast<int>(v.size());
    bool inside = false;
    for (int k = 0, j = n - 1; k < n; j = k++) {
      if (((v[k].second > cy) != (v[j].second > cy)) &&
          (cx <
           (v[j].first - v[k].first) * (cy - v[k].second) / (v[j].second - v[k].second + 1e-15) +
             v[k].first))
        inside = !inside;
    }
    if (inside)
      return true;
  }
  return false;
}

std::vector<Polymap::CDTEdge> Polymap::getCDTEdges(size_t max_edges) const {
  std::vector<CDTEdge> edges;
  if (!cdt_ready_ || max_edges == 0)
    return edges;

  // Do not reserve the full triangulation: this method is also used by the
  // ROS visualization path, where the caller may have a deliberately small
  // marker budget.
  const size_t reserve_hint = std::min<size_t>(max_edges, 1024);
  edges.reserve(reserve_hint);
  for (auto eit = cdt_.finite_edges_begin(); eit != cdt_.finite_edges_end(); ++eit) {
    if (edges.size() >= max_edges)
      break;
    auto face = eit->first;
    int idx = eit->second;
    auto pa = face->vertex(cdt_.cw(idx))->point();
    auto pb = face->vertex(cdt_.ccw(idx))->point();
    CDTEdge e;
    e.a = {static_cast<int>(pa.x()), static_cast<int>(pa.y())};
    e.b = {static_cast<int>(pb.x()), static_cast<int>(pb.y())};
    e.is_constrained = cdt_.is_constrained(*eit);
    edges.push_back(e);
  }
  return edges;
}

}  // namespace raystar
