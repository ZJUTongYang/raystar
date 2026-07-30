#pragma once

#include <cstddef>
#include <utility>

#include <raystar/cooperative_stop.h>

#include "visibility_validation.h"

namespace raystar {
namespace detail {

// Exact counterpart of the legacy pnpoly() return contract:
//   first  -> inside by the two-sided odd/even rule
//   second -> on the boundary, or a defensive left/right parity disagreement
//
// Visibility boundaries can contain same-ray far/near backtracking edges, so
// they are not assumed to satisfy CGAL::Polygon_2's simple-polygon
// precondition. Keeping both ray parities preserves the previous weakly-simple
// boundary semantics without division or double comparisons. Closed cycles
// use their implicit last->first edge. Open sectors instead use the real
// source->first and last->source spokes and never invent a last->first chord.
inline OperationStatus classifyPointInVisibilityBoundary(const VisibilityRegion& boundary,
                                                         const exact_geometry::Point& query,
                                                         const exact_geometry::Point& source,
                                                         VisibilityBoundaryMode mode,
                                                         std::pair<bool, bool>& classification,
                                                         const StopToken& stop_token) {
  if (stop_token.poll())
    return OperationStatus::stopped;

  if (mode == VisibilityBoundaryMode::closed_cycle && boundary.size() < 3)
    return OperationStatus::failure;
  if (mode == VisibilityBoundaryMode::open_sector && boundary.size() < 2)
    return OperationStatus::failure;
  if (mode != VisibilityBoundaryMode::closed_cycle && mode != VisibilityBoundaryMode::open_sector) {
    return OperationStatus::failure;
  }

  bool right_ray = false;
  bool left_ray = false;
  bool on_boundary = false;
  const auto process_edge = [&](const exact_geometry::Point& previous,
                                const exact_geometry::Point& current) {
    if (on_boundary)
      return;

    if (current == previous) {
      if (query == current)
        on_boundary = true;
      return;
    }

    if (exact_geometry::Kernel::Segment_2(previous, current).has_on(query)) {
      on_boundary = true;
      return;
    }

    const bool current_above = current.y() > query.y();
    const bool previous_above = previous.y() > query.y();
    if (current_above == previous_above)
      return;

    // Apply the same half-open y rule as the former implementation, but orient
    // every crossing from its lower endpoint to its upper endpoint. The y
    // denominator is then positive, so LEFT_TURN means the intersection lies
    // to the right of query and RIGHT_TURN means it lies to the left.
    const exact_geometry::Point& lower = current_above ? previous : current;
    const exact_geometry::Point& upper = current_above ? current : previous;
    const CGAL::Orientation side = CGAL::orientation(lower, upper, query);
    if (side == CGAL::LEFT_TURN)
      right_ray = !right_ray;
    else if (side == CGAL::RIGHT_TURN)
      left_ray = !left_ray;
    else
      on_boundary = true;
  };

  if (mode == VisibilityBoundaryMode::closed_cycle) {
    for (size_t i = 0, j = boundary.size() - 1; i < boundary.size(); j = i++) {
      if (stop_token.poll())
        return OperationStatus::stopped;
      process_edge(exactPoint(boundary[j]), exactPoint(boundary[i]));
    }
  } else {
    if (stop_token.poll())
      return OperationStatus::stopped;
    process_edge(source, exactPoint(boundary.front()));
    for (size_t i = 1; i < boundary.size(); ++i) {
      if (stop_token.poll())
        return OperationStatus::stopped;
      process_edge(exactPoint(boundary[i - 1]), exactPoint(boundary[i]));
    }
    if (stop_token.poll())
      return OperationStatus::stopped;
    process_edge(exactPoint(boundary.back()), source);
  }

  if (stop_token.poll())
    return OperationStatus::stopped;

  if (on_boundary)
    classification = {false, true};
  else if (right_ray == left_ray)
    classification = {right_ray, false};
  else
    classification = {false, true};
  return OperationStatus::success;
}

// Compatibility wrapper preserving the legacy result contract. Invalid input
// maps to {false, false}, as it did before cooperative stopping was added.
inline std::pair<bool, bool> classifyPointInVisibilityBoundary(const VisibilityRegion& boundary,
                                                               const exact_geometry::Point& query,
                                                               const exact_geometry::Point& source,
                                                               VisibilityBoundaryMode mode) {
  std::pair<bool, bool> classification{false, false};
  const StopToken never_stop;
  (void)classifyPointInVisibilityBoundary(
    boundary, query, source, mode, classification, never_stop);
  return classification;
}

// Cooperative overload for callers that supply an ordinary closed cycle.
inline OperationStatus classifyPointInVisibilityBoundary(const VisibilityRegion& boundary,
                                                         const exact_geometry::Point& query,
                                                         std::pair<bool, bool>& classification,
                                                         const StopToken& stop_token) {
  return classifyPointInVisibilityBoundary(
    boundary, query, query, VisibilityBoundaryMode::closed_cycle, classification, stop_token);
}

// Compatibility overload for callers that supply an ordinary closed cycle.
inline std::pair<bool, bool> classifyPointInVisibilityBoundary(const VisibilityRegion& boundary,
                                                               const exact_geometry::Point& query) {
  return classifyPointInVisibilityBoundary(
    boundary, query, query, VisibilityBoundaryMode::closed_cycle);
}

}  // namespace detail
}  // namespace raystar
