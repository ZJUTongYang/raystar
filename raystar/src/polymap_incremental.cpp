// Incremental occupancy updates for Polymap (R1 environment layer).
//
// Entry point: Polymap::applyOccupancyDelta.  The update applies newly
// occupied cells on top of an existing Polymap and rebuilds only what the
// change touches:
//
//   * Contours whose raw (pre-simplification) ring is unchanged keep their
//     obstacle indices and simplified geometry bit-for-bit ("frozen").
//   * Changed/merged contours are re-simplified under fresh appended
//     indices; the old indices become tombstones (empty vertex rings).
//   * A frozen contour whose simplified polygon is geometrically violated
//     by new geometry is "unfrozen": its raw ring is re-simplified under a
//     fresh index too (the unfreeze cascade).  This closes the gap that
//     simplification is conservative (polygon ⊇ cells), so new raw cells
//     can intrude into a frozen obstacle's historical fill area without
//     any grid adjacency.
//   * The CDT, vertex registry, and all validation gates are rebuilt in
//     full; any incremental-stage failure falls back to a full rebuild.
//
// Contract (see PolymapUpdateResult): old (obstacle, vertex) indices that
// are not retired keep pointing at bit-identical contours.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <raystar/polymap.h>

#include "polymap_detail.h"

namespace raystar {

namespace {

using polymap_impl::classifyExactSegments;
using polymap_impl::ExactSegmentRelation;
using polymap_impl::hasClockwiseOuterContour;
using polymap_impl::pointInPolygon;
using polymap_impl::validateCDTWithCGAL;

using Ring = std::vector<std::pair<int, int>>;

// Lexicographically smallest rotation of the ring.  Two extractions of the
// same obstacle produce the same directed cycle, but the starting vertex
// follows the extraction edge set order; normalizing the start makes ring
// equality well defined.  Direction is preserved: a reversed traversal of
// the same geometry simply fails to match and is re-simplified, which is
// fail-conservative, never wrong.
Ring canonicalRing(const Ring& ring) {
  if (ring.empty())
    return ring;
  const size_t size = ring.size();
  size_t best = 0;
  for (size_t start = 1; start < size; ++start) {
    for (size_t offset = 0; offset < size; ++offset) {
      const auto& best_point = ring[(best + offset) % size];
      const auto& start_point = ring[(start + offset) % size];
      if (start_point == best_point)
        continue;
      if (start_point < best_point)
        best = start;
      break;
    }
  }
  Ring canonical;
  canonical.reserve(size);
  for (size_t offset = 0; offset < size; ++offset)
    canonical.push_back(ring[(best + offset) % size]);
  return canonical;
}

std::uint64_t hashRing(const Ring& canonical) {
  std::uint64_t digest = 1469598103934665603ULL;
  for (const auto& point : canonical) {
    for (const int coordinate : {point.first, point.second}) {
      digest ^= static_cast<std::uint64_t>(static_cast<std::int32_t>(coordinate));
      digest *= 1099511628211ULL;
    }
  }
  return digest;
}

// Conservative geometric-conflict test between two simplified contours:
// any non-disjoint edge pair, or any vertex of one strictly inside the
// other, means the two obstacle polygons intrude into each other and one
// of them must be rebuilt.  Shared vertices between distinct obstacles
// are illegal in the registry, so any contact counts as a conflict.
//
// The outer contour is exempt from the containment half: its polygon is
// the free-space complement boundary (ray-cast "inside" covers the whole
// map), so containment against it is structural, not an intrusion.  Edge
// crossing remains meaningful for every pair including the outer contour.
struct ContourBox {
  int min_x = 0;
  int min_y = 0;
  int max_x = 0;
  int max_y = 0;
};

ContourBox contourBox(const Ring& ring) {
  ContourBox box{ring.front().first, ring.front().second, ring.front().first,
                 ring.front().second};
  for (const auto& point : ring) {
    box.min_x = std::min(box.min_x, point.first);
    box.max_x = std::max(box.max_x, point.first);
    box.min_y = std::min(box.min_y, point.second);
    box.max_y = std::max(box.max_y, point.second);
  }
  return box;
}

bool boxesOverlap(const ContourBox& first, const ContourBox& second) {
  return first.min_x <= second.max_x && second.min_x <= first.max_x &&
         first.min_y <= second.max_y && second.min_y <= first.max_y;
}

bool isClockwiseRing(const Ring& ring) {
  if (ring.size() < 3)
    return false;
  exact_geometry::FT twice_area(0);
  for (size_t index = 0; index < ring.size(); ++index) {
    const auto& from = ring[index];
    const auto& to = ring[(index + 1) % ring.size()];
    twice_area +=
      exact_geometry::FT(from.first) * to.second - exact_geometry::FT(to.first) * from.second;
  }
  return twice_area < exact_geometry::FT(0);
}

bool contoursConflict(const Ring& first,
                      const Ring& second,
                      bool either_is_outer,
                      const StopToken& stop_token) {
  if (first.size() < 2 || second.size() < 2)
    return false;
  if (!boxesOverlap(contourBox(first), contourBox(second)))
    return false;

  for (size_t first_edge = 0; first_edge < first.size(); ++first_edge) {
    if (stop_token.poll())
      return false;
    const auto& first_from = first[first_edge];
    const auto& first_to = first[(first_edge + 1) % first.size()];
    for (size_t second_edge = 0; second_edge < second.size(); ++second_edge) {
      const auto& second_from = second[second_edge];
      const auto& second_to = second[(second_edge + 1) % second.size()];
      const auto relation = classifyExactSegments(
        exact_geometry::Point(first_from.first, first_from.second),
        exact_geometry::Point(first_to.first, first_to.second),
        exact_geometry::Point(second_from.first, second_from.second),
        exact_geometry::Point(second_to.first, second_to.second));
      if (relation.relation != ExactSegmentRelation::disjoint)
        return true;
    }
  }

  if (either_is_outer)
    return stop_token.poll();

  // Containment without edge crossing: a new contour can sit entirely (or
  // with all vertices) inside a frozen obstacle's historical fill area.
  // Ray-cast containment over conservatively-simplified rings may also
  // flag a legitimately encircled obstacle (an island inside a U-shaped
  // wall's embrace); that over-unfreezes but never corrupts -- the re-
  // simplification reproduces the same legal contour.
  const auto anyVertexInside = [&](const Ring& outer, const Ring& points) {
    for (const auto& point : points) {
      bool inside = false;
      if (pointInPolygon(outer, point.first, point.second, stop_token, inside) !=
          OperationStatus::success)
        return false;  // stopped; the caller re-checks the stop token
      if (inside)
        return true;
    }
    return false;
  };
  if (anyVertexInside(first, second) || anyVertexInside(second, first))
    return true;
  return stop_token.poll();
}

// Index of the first clockwise (negative signed area) contour, mirroring
// hasClockwiseOuterContour's outer-contour detection.  Empty (tombstone)
// rings have zero area and are skipped.
int clockwiseContourIndex(const std::vector<Obs>& obstacles) {
  for (size_t index = 0; index < obstacles.size(); ++index) {
    const auto& ring = obstacles[index].ordered_vertices_;
    if (ring.size() < 3)
      continue;
    exact_geometry::FT twice_area(0);
    for (size_t vertex = 0; vertex < ring.size(); ++vertex) {
      const auto& from = ring[vertex];
      const auto& to = ring[(vertex + 1) % ring.size()];
      twice_area +=
        exact_geometry::FT(from.first) * to.second - exact_geometry::FT(to.first) * from.second;
    }
    if (twice_area < exact_geometry::FT(0))
      return static_cast<int>(index);
  }
  return -1;
}

}  // namespace

PolymapUpdateResult Polymap::applyOccupancyDelta(const Polymap& base,
                                                 const std::vector<std::pair<int, int>>& newly_occupied_cells,
                                                 int start_x,
                                                 int start_y,
                                                 const Point2d& start_position,
                                                 const std::vector<PolymapEndpoint>& goals,
                                                 const StopToken& stop_token,
                                                 const PlanningLimits& limits) {
  PolymapUpdateResult result;
  if (stop_token.poll()) {
    result.status = PolymapCreateStatus::stopped;
    return result;
  }
  if (base.xsize_ <= 0 || base.ysize_ <= 0) {
    result.error = "Base polymap has no occupancy";
    return result;
  }

  // Snapshot the post-update occupancy once: the incremental path works on
  // it and any fallback rebuilds from it, so both see the same world.
  std::vector<uint8_t> updated = base.data_;
  bool changed = false;
  for (const auto& cell : newly_occupied_cells) {
    if (cell.first < 0 || cell.second < 0 || cell.first >= base.xsize_ ||
        cell.second >= base.ysize_) {
      result.error = "Newly occupied cell (" + std::to_string(cell.first) + ", " +
                     std::to_string(cell.second) + ") is outside the base occupancy";
      return result;
    }
    const size_t index =
      static_cast<size_t>(cell.second) * static_cast<size_t>(base.xsize_) +
      static_cast<size_t>(cell.first);
    if (updated[index] == 0) {
      updated[index] = 1;
      changed = true;
    }
  }
  if (!changed) {
    result.error = "No newly occupied cell changes the base occupancy";
    return result;
  }

  const auto fall_back = [&](const std::string& reason) {
    GridMap rebuilt;
    rebuilt.width = static_cast<unsigned int>(base.xsize_);
    rebuilt.height = static_cast<unsigned int>(base.ysize_);
    rebuilt.resolution = 1.0f;
    rebuilt.origin_x = 0.0;
    rebuilt.origin_y = 0.0;
    rebuilt.data = updated;
    auto rebuild =
      Polymap::create(rebuilt, start_x, start_y, start_position, goals, stop_token, limits);
    result.status = rebuild.status;
    result.error = rebuild.error;
    result.fallback_reason = reason;
    // Reset the partial change lists from the declined incremental path:
    // on a successful rebuild they cover every live obstacle on both
    // sides; on a failed one they must not linger half-filled.
    result.retired_obstacles.clear();
    result.added_obstacles.clear();
    if (rebuild) {
      result.value = std::move(rebuild.value);
      result.fell_back_to_full_rebuild = true;
      for (size_t index = 0; index < base.obs_.size(); ++index) {
        if (!base.obs_[index].ordered_vertices_.empty())
          result.retired_obstacles.push_back(static_cast<int>(index));
      }
      if (result.value) {
        const auto& obstacles = result.value->obstacles();
        for (size_t index = 0; index < obstacles.size(); ++index) {
          if (!obstacles[index].ordered_vertices_.empty())
            result.added_obstacles.push_back(static_cast<int>(index));
        }
      }
    }
  };

  bool declined = false;
  Polymap candidate(base,
                    newly_occupied_cells,
                    start_x,
                    start_y,
                    start_position,
                    goals,
                    stop_token,
                    result.retired_obstacles,
                    result.added_obstacles,
                    declined);
  if (stop_token.poll() && candidate.construction_stopped_) {
    result.status = PolymapCreateStatus::stopped;
    return result;
  }
  if (declined) {
    // The incremental path reports why it declined through the candidate's
    // construction error only for diagnostics; the fallback decides the
    // final outcome.
    fall_back(candidate.construction_error_.empty()
                ? std::string("incremental update declined")
                : candidate.construction_error_);
    return result;
  }
  if (candidate.no_path_) {
    // Same shape as Polymap::create on no_path: status only, no value.
    result.status = PolymapCreateStatus::no_path;
    return result;
  }
  if (!candidate.solution_exist_ || !candidate.cdt_ready_) {
    result.error = candidate.construction_error_.empty()
                     ? std::string("Incremental polymap construction failed")
                     : candidate.construction_error_;
    return result;
  }
  result.status = PolymapCreateStatus::ready;
  result.value = std::move(candidate);
  return result;
}

Polymap::Polymap(const Polymap& base,
                 const std::vector<std::pair<int, int>>& newly_occupied_cells,
                 int start_x,
                 int start_y,
                 const Point2d& start_position,
                 const std::vector<PolymapEndpoint>& goals,
                 const StopToken& stop_token,
                 std::vector<int>& retired,
                 std::vector<int>& added,
                 bool& declined)
  : xsize_(base.xsize_), ysize_(base.ysize_) {
  const auto decline = [&](const std::string& message) {
    declined = true;
    construction_error_ = message;
    clearStoppedConstructionState();
  };

  if (stop_token.poll()) {
    clearStoppedConstructionState();
    return;
  }

  // ---- Stage 1: occupancy + endpoint admission (same rules as create) --
  const size_t cell_count = static_cast<size_t>(xsize_) * static_cast<size_t>(ysize_);
  data_ = base.data_;
  for (const auto& cell : newly_occupied_cells) {
    const size_t index =
      static_cast<size_t>(cell.second) * static_cast<size_t>(xsize_) +
      static_cast<size_t>(cell.first);
    data_[index] = 1;
  }
  vertices_location_x_flat_.resize(cell_count, -1);
  vertices_location_y_flat_.resize(cell_count, -1);

  if (!std::isfinite(start_position.first) || !std::isfinite(start_position.second)) {
    construction_error_ = "Start position is not finite";
    return;
  }
  if (goals.empty()) {
    construction_error_ = "At least one goal is required";
    return;
  }
  if (std::floor(start_position.first) != static_cast<double>(start_x) ||
      std::floor(start_position.second) != static_cast<double>(start_y)) {
    construction_error_ = "Start continuous position does not belong to the supplied grid cell";
    return;
  }
  const size_t start_index =
    static_cast<size_t>(start_x) + static_cast<size_t>(start_y) * static_cast<size_t>(xsize_);
  if (start_index >= data_.size() || data_[start_index] != 0) {
    construction_error_ = "Start grid cell is occupied";
    return;
  }
  for (const auto& goal : goals) {
    if (std::floor(goal.position.first) != static_cast<double>(goal.cell_x) ||
        std::floor(goal.position.second) != static_cast<double>(goal.cell_y)) {
      construction_error_ = "Goal continuous position does not belong to the supplied grid cell";
      return;
    }
    const size_t goal_index = static_cast<size_t>(goal.cell_x) +
                              static_cast<size_t>(goal.cell_y) * static_cast<size_t>(xsize_);
    if (goal_index >= data_.size() || data_[goal_index] != 0) {
      construction_error_ = "Goal grid cell is occupied";
      return;
    }
  }

  // ---- Stage 2: full re-extraction of raw contours (cheap, O(W*H)) -----
  // getPolyObstacles fills obs_ with the new raw rings and captures
  // raw_obstacles_ alongside them.
  const OperationStatus obstacle_status = getPolyObstacles(start_x, start_y, goals, stop_token);
  if (obstacle_status == OperationStatus::stopped) {
    clearStoppedConstructionState();
    return;
  }
  solution_exist_ = obstacle_status == OperationStatus::success;
  if (!solution_exist_) {
    no_path_ = true;
    construction_error_ = "Start and every goal must be in the same reachable free-space component";
    return;
  }

  // ---- Stage 3: ring matching and tombstone assembly -------------------
  // Match each new raw ring against the base's raw rings; matched rings
  // freeze (reuse the base's simplified contour under the same index),
  // unmatched ones are appended for simplification, and base rings without
  // a match become tombstones.
  std::vector<Obs> new_rings = std::move(obs_);
  std::vector<Ring> new_raw = std::move(raw_obstacles_);
  obs_.clear();
  raw_obstacles_.clear();

  std::unordered_map<std::uint64_t, std::vector<int>> base_by_hash;
  // canonical ring cached per base index: the bucket is built once and the
  // equality check below reuses it instead of recomputing the rotation.
  std::vector<Ring> base_canonical(base.raw_obstacles_.size());
  for (size_t index = 0; index < base.raw_obstacles_.size(); ++index) {
    const auto& ring = base.raw_obstacles_[index];
    if (ring.empty())
      continue;  // historical tombstone stays a tombstone
    base_canonical[index] = canonicalRing(ring);
    base_by_hash[hashRing(base_canonical[index])].push_back(static_cast<int>(index));
  }

  const size_t base_obstacle_count = base.obs_.size();
  std::vector<int> match_of_new(new_rings.size(), -1);
  std::vector<char> base_matched(base_obstacle_count, 0);
  for (size_t index = 0; index < new_rings.size(); ++index) {
    if (stop_token.poll()) {
      clearStoppedConstructionState();
      return;
    }
    const Ring canonical = canonicalRing(new_raw[index]);
    const auto found = base_by_hash.find(hashRing(canonical));
    if (found == base_by_hash.end())
      continue;
    for (const int candidate_index : found->second) {
      if (base_matched[candidate_index])
        continue;
      if (base_canonical[candidate_index] == canonical) {
        match_of_new[index] = candidate_index;
        base_matched[candidate_index] = 1;
        break;
      }
    }
  }

  size_t appended = 0;
  for (size_t index = 0; index < new_rings.size(); ++index) {
    if (match_of_new[index] < 0)
      ++appended;
  }

  std::vector<Obs> assembled(base_obstacle_count + appended);
  std::vector<Ring> assembled_raw(base_obstacle_count + appended);
  for (size_t index = 0; index < base_obstacle_count; ++index) {
    if (base_matched[index]) {
      assembled[index] = base.obs_[index];
      assembled_raw[index] = base.raw_obstacles_[index];
    } else if (!base.raw_obstacles_[index].empty()) {
      retired.push_back(static_cast<int>(index));  // tombstone: keep slot, empty ring
    }
  }
  std::vector<size_t> to_simplify;
  size_t next_append = base_obstacle_count;
  for (size_t index = 0; index < new_rings.size(); ++index) {
    if (match_of_new[index] >= 0)
      continue;
    assembled[next_append].ordered_vertices_ = std::move(new_rings[index].ordered_vertices_);
    assembled_raw[next_append] = std::move(new_raw[index]);
    added.push_back(static_cast<int>(next_append));
    to_simplify.push_back(next_append);
    ++next_append;
  }
  obs_ = std::move(assembled);
  raw_obstacles_ = std::move(assembled_raw);

  // ---- Stage 4: outer-contour protection (V1: refuse and let the -------
  // caller fall back) ----------------------------------------------------
  const int base_outer = clockwiseContourIndex(base.obs_);
  if (base_outer >= 0 && !base_matched[base_outer]) {
    decline("Incremental update would affect the outer contour");
    return;
  }
  const bool has_outer_contour = hasClockwiseOuterContour(obs_, stop_token);
  if (stop_token.poll()) {
    clearStoppedConstructionState();
    return;
  }

  // ---- Stage 5: endpoint interiority against the pre-simplify rings ----
  std::string endpoint_error;
  if (has_outer_contour) {
    OperationStatus status = validateFreeSpaceInterior(start_position, stop_token, &endpoint_error);
    if (status == OperationStatus::stopped) {
      clearStoppedConstructionState();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid start position: " + endpoint_error;
      return;
    }
    for (const auto& goal : goals) {
      status = validateFreeSpaceInterior(goal.position, stop_token, &endpoint_error);
      if (status == OperationStatus::stopped) {
        clearStoppedConstructionState();
        return;
      }
      if (status == OperationStatus::failure) {
        construction_error_ = "Invalid goal position: " + endpoint_error;
        return;
      }
    }
  }
  OperationStatus status = validateObstacleTopology(construction_error_, false, stop_token);
  if (status == OperationStatus::stopped) {
    clearStoppedConstructionState();
    return;
  }
  if (status == OperationStatus::failure) {
    declined = true;  // malformed extraction: rebuild from scratch
    return;
  }

  // ---- Stage 6: simplify only appended contours -------------------------
  std::vector<Point2d> protected_points;
  protected_points.reserve(goals.size() + 1);
  protected_points.emplace_back(start_position);
  for (const auto& goal : goals)
    protected_points.emplace_back(goal.position);
  if (!to_simplify.empty()) {
    if (!simplifyPolyObstaclesImpl(protected_points, stop_token, &to_simplify)) {
      if (stop_token.poll()) {
        clearStoppedConstructionState();
        return;
      }
      declined = true;
      construction_error_ = "Incremental obstacle simplification failed";
      return;
    }
  }

  // ---- Stage 7: unfreeze cascade ----------------------------------------
  // Frozen contours keep geometry that was simplified under the old global
  // context; new geometry may intrude into their conservative fill area.
  // Any frozen contour conflicting with appended geometry is unfrozen
  // (tombstoned and re-simplified under a fresh index) until no conflicts
  // remain.  Declining is budgeted: if nearly everything unfreezes, a full
  // rebuild is the honest result.
  size_t live_frozen = 0;
  for (size_t index = 0; index < base_obstacle_count; ++index) {
    if (base_matched[index])
      ++live_frozen;
  }
  // A proportional budget: when a change would unfreeze more than half of
  // the frozen contours, a full rebuild is the honest result (and avoids
  // the O(victims x frozen x appended) rescan of a near-total cascade).
  const size_t unfreeze_budget = std::max<size_t>(1, live_frozen / 2);
  size_t unfrozen = 0;
  std::vector<size_t> simplify_batch;
  while (true) {
    if (stop_token.poll()) {
      clearStoppedConstructionState();
      return;
    }
    int victim = -1;
    for (size_t frozen_index = 0; frozen_index < base_obstacle_count; ++frozen_index) {
      if (!base_matched[frozen_index])
        continue;
      const auto& frozen_ring = obs_[frozen_index].ordered_vertices_;
      if (frozen_ring.empty())
        continue;  // unfrozen earlier in this cascade
      bool conflict = false;
      for (const int appended_index : added) {
        const auto& appended_ring =
          obs_[static_cast<size_t>(appended_index)].ordered_vertices_;
        if (appended_ring.empty())
          continue;
        const bool either_is_outer =
          isClockwiseRing(frozen_ring) || isClockwiseRing(appended_ring);
        if (contoursConflict(frozen_ring, appended_ring, either_is_outer, stop_token)) {
          conflict = true;
          break;
        }
      }
      if (stop_token.poll()) {
        clearStoppedConstructionState();
        return;
      }
      if (conflict) {
        victim = static_cast<int>(frozen_index);
        break;
      }
    }
    if (victim < 0)
      break;
    if (++unfrozen >= unfreeze_budget) {
      decline("Incremental unfreeze cascade exceeded its budget");
      return;
    }
    // Unfreeze: tombstone the old slot, append a re-simplification of its
    // raw ring under a fresh index.
    const Ring unfrozen_raw = raw_obstacles_[static_cast<size_t>(victim)];
    retired.push_back(victim);
    obs_[static_cast<size_t>(victim)].ordered_vertices_.clear();
    raw_obstacles_[static_cast<size_t>(victim)].clear();
    base_matched[static_cast<size_t>(victim)] = 0;
    const size_t fresh = obs_.size();
    Obs fresh_obstacle;
    fresh_obstacle.ordered_vertices_ = unfrozen_raw;
    obs_.push_back(std::move(fresh_obstacle));
    raw_obstacles_.push_back(unfrozen_raw);
    added.push_back(static_cast<int>(fresh));
    simplify_batch.push_back(fresh);
  }
  if (!simplify_batch.empty()) {
    if (!simplifyPolyObstaclesImpl(protected_points, stop_token, &simplify_batch)) {
      if (stop_token.poll()) {
        clearStoppedConstructionState();
        return;
      }
      decline("Unfreeze-cascade simplification failed");
      return;
    }
  }

  // ---- Stage 8: validation gates and full downstream rebuild ------------
  status = validateObstacleTopology(construction_error_, true, stop_token);
  if (status == OperationStatus::stopped) {
    clearStoppedConstructionState();
    return;
  }
  if (status == OperationStatus::failure) {
    decline("Incremental update rejected by the edge-relation validator");
    return;
  }
  if (has_outer_contour) {
    status = validateFreeSpaceInterior(start_position, stop_token, &endpoint_error);
    if (status == OperationStatus::stopped) {
      clearStoppedConstructionState();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid start position: " + endpoint_error;
      return;
    }
    for (const auto& goal : goals) {
      status = validateFreeSpaceInterior(goal.position, stop_token, &endpoint_error);
      if (status == OperationStatus::stopped) {
        clearStoppedConstructionState();
        return;
      }
      if (status == OperationStatus::failure) {
        construction_error_ = "Invalid goal position: " + endpoint_error;
        return;
      }
    }
  }
  status = registerVertices(construction_error_, stop_token);
  if (status == OperationStatus::stopped) {
    clearStoppedConstructionState();
    return;
  }
  if (status == OperationStatus::failure) {
    decline("Incremental update rejected by the vertex registry");
    return;
  }
  status = constructCGALRelated(&validateCDTWithCGAL, construction_error_, stop_token);
  if (status == OperationStatus::stopped) {
    clearStoppedConstructionState();
    return;
  }
  if (status == OperationStatus::failure) {
    decline("Incremental update rejected by the CDT construction");
    return;
  }
}

}  // namespace raystar
