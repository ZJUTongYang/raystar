#include "visibility_validation.h"

#include <cmath>
#include <string>
#include <utility>
#include <vector>

namespace raystar {
namespace detail {
namespace {

using ExactPoint = exact_geometry::Point;

struct ExactEdge {
  ExactPoint from;
  ExactPoint to;
};

bool areStrictlyOpposite(CGAL::Orientation first, CGAL::Orientation second) {
  return (first == CGAL::LEFT_TURN && second == CGAL::RIGHT_TURN) ||
         (first == CGAL::RIGHT_TURN && second == CGAL::LEFT_TURN);
}

bool hasProperInteriorCrossing(const ExactEdge& first, const ExactEdge& second) {
  if (first.from == first.to || second.from == second.to)
    return false;

  const CGAL::Orientation first_from = CGAL::orientation(first.from, first.to, second.from);
  const CGAL::Orientation first_to = CGAL::orientation(first.from, first.to, second.to);
  const CGAL::Orientation second_from = CGAL::orientation(second.from, second.to, first.from);
  const CGAL::Orientation second_to = CGAL::orientation(second.from, second.to, first.to);

  return areStrictlyOpposite(first_from, first_to) && areStrictlyOpposite(second_from, second_to);
}

bool areAdjacentCycleEdges(size_t first, size_t second, size_t edge_count) {
  return (first + 1) % edge_count == second || (second + 1) % edge_count == first;
}

}  // namespace

bool validateVisibilityGeometry(const VisibilityRegion& visibility_region,
                                const VisibilityGeometryContext& context,
                                std::string* error) {
  const StopToken never_stop;
  return validateVisibilityGeometry(visibility_region, context, never_stop, error) ==
         OperationStatus::success;
}

OperationStatus validateVisibilityGeometry(const VisibilityRegion& visibility_region,
                                           const VisibilityGeometryContext& context,
                                           const StopToken& stop_token,
                                           std::string* error) {
  if (error)
    error->clear();
  const auto fail = [&](std::string message) {
    if (error)
      *error = std::move(message);
    return OperationStatus::failure;
  };
  const auto stopRequested = [&]() {
    if (!stop_token.poll())
      return false;
    if (error)
      error->clear();
    return true;
  };

  if (stopRequested())
    return OperationStatus::stopped;

  if (context.mode != VisibilityBoundaryMode::closed_cycle &&
      context.mode != VisibilityBoundaryMode::open_sector) {
    return fail("Visibility boundary has an unsupported closure mode");
  }
  if (visibility_region.size() < 2)
    return fail("Visibility boundary has fewer than two endpoints");
  if (context.mode == VisibilityBoundaryMode::closed_cycle && visibility_region.size() < 3) {
    return fail("Closed visibility boundary has fewer than three endpoints");
  }

  if (context.mode == VisibilityBoundaryMode::open_sector) {
    if (!context.start_anchor || !context.end_anchor)
      return fail("Open visibility boundary is missing a sector anchor");
    if (*context.start_anchor == context.source || *context.end_anchor == context.source) {
      return fail("Open visibility sector anchor coincides with source");
    }
  } else if (context.start_anchor || context.end_anchor) {
    return fail("Closed visibility boundary unexpectedly has a sector anchor");
  }

  std::vector<ExactPoint> points;
  points.reserve(visibility_region.size());
  for (size_t index = 0; index < visibility_region.size(); ++index) {
    if (stopRequested())
      return OperationStatus::stopped;

    const auto& endpoint = visibility_region[index];
    if (!std::isfinite(endpoint.position.first) || !std::isfinite(endpoint.position.second)) {
      return fail("Visibility endpoint " + std::to_string(index) + " is not finite");
    }

    points.emplace_back(exactPoint(endpoint));
    if (points.back() == context.source) {
      return fail("Visibility endpoint " + std::to_string(index) + " coincides with source");
    }
  }

  if (context.mode == VisibilityBoundaryMode::open_sector &&
      (points.front() != *context.start_anchor || points.back() != *context.end_anchor)) {
    return fail("Open visibility boundary does not match its sector anchors");
  }

  // An explicit visibility edge may touch or retrace another edge, but it may
  // not run through the source.  Zero-length edges are deliberately accepted:
  // existing consumers treat an exact repeated endpoint as an ignorable
  // discontinuity.
  const size_t explicit_edge_count =
    context.mode == VisibilityBoundaryMode::closed_cycle ? points.size() : points.size() - 1;
  for (size_t index = 0; index < explicit_edge_count; ++index) {
    if (stopRequested())
      return OperationStatus::stopped;

    const ExactPoint& from = points[index];
    const ExactPoint& to = context.mode == VisibilityBoundaryMode::closed_cycle
                             ? points[(index + 1) % points.size()]
                             : points[index + 1];
    if (from != to && exact_geometry::Kernel::Segment_2(from, to).has_on(context.source)) {
      return fail("Visibility boundary edge " + std::to_string(index) + " passes through source");
    }
  }

  // Use the boundary's real closure semantics for crossing checks.  An open
  // sector is closed conceptually as source -> first -> ... -> last -> source;
  // the artificial last -> first chord is not part of its boundary.
  std::vector<ExactPoint> cycle_vertices;
  cycle_vertices.reserve(points.size() + 1);
  if (context.mode == VisibilityBoundaryMode::open_sector)
    cycle_vertices.emplace_back(context.source);
  for (const auto& point : points) {
    if (stopRequested())
      return OperationStatus::stopped;
    cycle_vertices.emplace_back(point);
  }

  std::vector<ExactEdge> cycle_edges;
  cycle_edges.reserve(cycle_vertices.size());
  for (size_t index = 0; index < cycle_vertices.size(); ++index) {
    if (stopRequested())
      return OperationStatus::stopped;

    cycle_edges.push_back(
      ExactEdge{cycle_vertices[index], cycle_vertices[(index + 1) % cycle_vertices.size()]});
  }

  size_t crossing_pair_count = 0;
  for (size_t first = 0; first < cycle_edges.size(); ++first) {
    for (size_t second = first + 1; second < cycle_edges.size(); ++second) {
      if (crossing_pair_count++ % 64 == 0 && stopRequested())
        return OperationStatus::stopped;

      if (areAdjacentCycleEdges(first, second, cycle_edges.size()))
        continue;
      if (hasProperInteriorCrossing(cycle_edges[first], cycle_edges[second])) {
        return fail("Visibility boundary has a proper crossing between edges " +
                    std::to_string(first) + " and " + std::to_string(second));
      }
    }
  }

  if (context.mode == VisibilityBoundaryMode::closed_cycle) {
    for (size_t index = 1; index < points.size(); ++index) {
      if (stopRequested())
        return OperationStatus::stopped;

      if (exact_geometry::compareRaysLikeAtan2(context.source, points[index - 1], points[index]) ==
          CGAL::LARGER) {
        return fail("Visibility ray order decreases between endpoints " +
                    std::to_string(index - 1) + " and " + std::to_string(index));
      }
    }
    if (stopRequested())
      return OperationStatus::stopped;
    return OperationStatus::success;
  }

  const ExactPoint& start = *context.start_anchor;
  const ExactPoint& end = *context.end_anchor;
  for (size_t index = 0; index < points.size(); ++index) {
    if (stopRequested())
      return OperationStatus::stopped;

    if (!exact_geometry::isClosedCounterClockwiseSweepMember(
          context.source, start, end, points[index])) {
      return fail("Visibility endpoint " + std::to_string(index) +
                  " lies outside the requested sector");
    }
    if (index > 0 && exact_geometry::compareRaysCounterClockwiseFrom(
                       context.source, start, points[index - 1], points[index]) == CGAL::LARGER) {
      return fail("Visibility ray order decreases between endpoints " + std::to_string(index - 1) +
                  " and " + std::to_string(index));
    }
  }
  if (stopRequested())
    return OperationStatus::stopped;
  return OperationStatus::success;
}

}  // namespace detail
}  // namespace raystar
