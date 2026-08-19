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

#include "polymap_detail.h"

namespace raystar {

using namespace polymap_impl;

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

}  // namespace raystar
