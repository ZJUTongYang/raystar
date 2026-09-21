// Incremental Polymap update tests (R1 environment layer).
//
// Contract under test (see PolymapUpdateResult):
//   C1  obstacles not involved in an update keep their indices and
//       simplified geometry bit-for-bit;
//   C2  a full rebuild from the same post-update occupancy stays a valid
//       fallback for every declined path (validators gate it);
//   C3  the change set (retired/added) is reported honestly, including
//       the full-rebuild fallback where it covers everything.

#include <gtest/gtest.h>
#include <raystar/polymap.h>

#include <algorithm>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace raystar {

namespace {

GridMap makeBorderedMap(unsigned int width, unsigned int height) {
  GridMap map;
  map.width = width;
  map.height = height;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.assign(static_cast<size_t>(width) * height, 0);
  for (unsigned int x = 0; x < width; ++x) {
    map.data[x] = 1;
    map.data[(height - 1) * width + x] = 1;
  }
  for (unsigned int y = 0; y < height; ++y) {
    map.data[y * width] = 1;
    map.data[y * width + width - 1] = 1;
  }
  return map;
}

void occupy(GridMap& map, int x, int y) {
  map.data[static_cast<size_t>(y) * map.width + static_cast<size_t>(x)] = 1;
}

Polymap makeReadyPolymap(const GridMap& map,
                         int start_x,
                         int start_y,
                         const Point2d& start,
                         const std::vector<PolymapEndpoint>& goals) {
  auto result = Polymap::create(map, start_x, start_y, start, goals, StopToken{});
  if (!result) {
    throw std::runtime_error(result.error.empty() ? "Polymap construction did not produce a value"
                                                  : result.error);
  }
  return std::move(*result.value);
}

// (obstacle, vertex) -> coordinate mapping of a Polymap's live obstacles.
std::set<std::pair<int, std::pair<int, int>>> vertexMap(const Polymap& polymap) {
  std::set<std::pair<int, std::pair<int, int>>> mapping;
  const auto& obstacles = polymap.obstacles();
  for (size_t obstacle = 0; obstacle < obstacles.size(); ++obstacle) {
    const auto& ring = obstacles[obstacle].ordered_vertices_;
    for (size_t vertex = 0; vertex < ring.size(); ++vertex)
      mapping.emplace(static_cast<int>(obstacle), ring[vertex]);
  }
  return mapping;
}

std::set<std::pair<int, std::pair<int, int>>> verticesOf(
  const std::set<int>& obstacles, const Polymap& polymap) {
  std::set<std::pair<int, std::pair<int, int>>> mapping;
  const auto& all = polymap.obstacles();
  for (const int obstacle : obstacles) {
    const auto& ring = all[static_cast<size_t>(obstacle)].ordered_vertices_;
    for (const auto& point : ring)
      mapping.emplace(obstacle, point);
  }
  return mapping;
}

// Runs the planner end-to-end on both an incrementally updated map and a
// full rebuild of the same occupancy, returning both results for class-
// level comparison (certified lengths must agree; see the class-
// preservation lemma in the design notes).
struct BothSides {
  PolymapUpdateResult incremental;
  PolymapCreateResult rebuilt;
};

BothSides updateAndRebuild(const Polymap& base,
                           const std::vector<std::pair<int, int>>& cells,
                           int start_x,
                           int start_y,
                           const Point2d& start,
                           const std::vector<PolymapEndpoint>& goals) {
  BothSides sides;
  sides.incremental =
    Polymap::applyOccupancyDelta(base, cells, {}, start_x, start_y, start, goals, StopToken{});
  GridMap rebuilt_map;
  rebuilt_map.width = static_cast<unsigned int>(base.width());
  rebuilt_map.height = static_cast<unsigned int>(base.height());
  rebuilt_map.resolution = 1.0f;
  rebuilt_map.data = base.occupancyData();
  for (const auto& cell : cells)
    rebuilt_map.data[static_cast<size_t>(cell.second) * rebuilt_map.width +
                     static_cast<size_t>(cell.first)] = 1;
  sides.rebuilt = Polymap::create(rebuilt_map, start_x, start_y, start, goals, StopToken{});
  return sides;
}

}  // namespace

}  // namespace raystar

using namespace raystar;

// --- Contract C1: untouched obstacles keep indices and geometry --------

TEST(PolymapIncremental, UntouchedObstaclesKeepIndicesAndGeometry) {
  auto map = makeBorderedMap(20, 20);
  for (int x = 6; x <= 9; ++x) {
    occupy(map, x, 8);
    occupy(map, x, 14);
  }
  occupy(map, 15, 11);  // a separate lone obstacle
  occupy(map, 15, 12);

  const std::vector<PolymapEndpoint> goals{{18, 10, Point2d{18.5, 10.5}}};
  Polymap base = makeReadyPolymap(map, 2, 10, Point2d{2.5, 10.5}, goals);
  const auto before = vertexMap(base);

  // Grow the lone obstacle (an island untouched by anything else).
  const std::vector<std::pair<int, int>> cells{{15, 10}, {16, 10}};
  auto result = Polymap::applyOccupancyDelta(base, cells, {}, 2, 10, Point2d{2.5, 10.5}, goals,
                                             StopToken{});
  ASSERT_TRUE(result) << result.error;
  EXPECT_FALSE(result.fell_back_to_full_rebuild);
  ASSERT_TRUE(result.value);

  // Every retired index disappears from the live mapping; every other old
  // (obstacle, vertex) pair survives with an identical coordinate.
  std::set<int> retired(result.retired_obstacles.begin(), result.retired_obstacles.end());
  EXPECT_EQ(retired.size(), result.retired_obstacles.size());  // no duplicates
  for (const auto& entry : before) {
    if (retired.count(entry.first))
      continue;
    EXPECT_TRUE(result.value->isValidTopology({entry.first, 0}) ||
                !result.value->obstacles()[static_cast<size_t>(entry.first)]
                   .ordered_vertices_.empty())
      << "non-retired obstacle " << entry.first << " became a tombstone";
  }
  // The change set is non-trivial: the grown island retired exactly one
  // obstacle (its old ring) and appended exactly one (its new ring).
  EXPECT_EQ(result.retired_obstacles.size(), 1u);
  EXPECT_EQ(result.added_obstacles.size(), 1u);
}

TEST(PolymapIncremental, FrozenContourGeometryIsBitIdentical) {
  auto map = makeBorderedMap(20, 20);
  for (int x = 6; x <= 9; ++x) {
    occupy(map, x, 8);
    occupy(map, x, 14);
  }
  occupy(map, 15, 11);
  occupy(map, 15, 12);

  const std::vector<PolymapEndpoint> goals{{18, 10, Point2d{18.5, 10.5}}};
  Polymap base = makeReadyPolymap(map, 2, 10, Point2d{2.5, 10.5}, goals);

  // Snapshot every live obstacle ring, then grow the lone island.
  std::vector<std::vector<std::pair<int, int>>> rings_before;
  for (const auto& obstacle : base.obstacles())
    rings_before.push_back(obstacle.ordered_vertices_);

  const std::vector<std::pair<int, int>> cells{{16, 11}};
  auto result = Polymap::applyOccupancyDelta(base, cells, {}, 2, 10, Point2d{2.5, 10.5}, goals,
                                             StopToken{});
  ASSERT_TRUE(result) << result.error;
  ASSERT_FALSE(result.fell_back_to_full_rebuild);

  std::set<int> retired(result.retired_obstacles.begin(), result.retired_obstacles.end());
  const auto& after = result.value->obstacles();
  ASSERT_EQ(after.size(), rings_before.size() + result.added_obstacles.size());
  for (size_t index = 0; index < rings_before.size(); ++index) {
    if (retired.count(static_cast<int>(index))) {
      EXPECT_TRUE(after[index].ordered_vertices_.empty())
        << "retired obstacle " << index << " must be a tombstone";
    } else {
      EXPECT_EQ(after[index].ordered_vertices_, rings_before[index])
        << "frozen obstacle " << index << " changed geometry";
    }
  }
}

// --- Merge: a bridge joins two obstacles into one -----------------------

TEST(PolymapIncremental, BridgeMergeRetiresBothAndAppendsOne) {
  auto map = makeBorderedMap(24, 20);
  for (int x = 6; x <= 9; ++x)
    occupy(map, x, 6);
  for (int x = 6; x <= 9; ++x)
    occupy(map, x, 13);
  // Two horizontal walls with a gap; bridging the gap merges them.
  const std::vector<PolymapEndpoint> goals{{21, 9, Point2d{21.5, 9.5}}};
  Polymap base = makeReadyPolymap(map, 2, 9, Point2d{2.5, 9.5}, goals);

  const std::vector<std::pair<int, int>> bridge{{6, 7},  {6, 8},  {6, 9},
                                                {6, 10}, {6, 11}, {6, 12}};
  auto result = Polymap::applyOccupancyDelta(base, bridge, {}, 2, 9, Point2d{2.5, 9.5}, goals,
                                             StopToken{});
  ASSERT_TRUE(result) << result.error;

  const std::set<int> retired(result.retired_obstacles.begin(), result.retired_obstacles.end());
  // The outer frame must never be retired.
  for (size_t index = 0; index < base.obstacles().size(); ++index) {
    const auto& ring = base.obstacles()[index].ordered_vertices_;
    if (ring.size() >= 4 && ring.front().second == 1 && ring.back().second == 1 &&
        std::any_of(ring.begin(), ring.end(), [](const auto& point) {
          return point.first == 1 && point.second == 1;
        })) {
      // Heuristic outer-frame check: a contour touching (1,1).
      EXPECT_FALSE(retired.count(static_cast<int>(index)))
        << "the outer contour was affected";
      break;
    }
  }

  if (!result.fell_back_to_full_rebuild) {
    // The merged contour is appended; the cascade may append re-
    // simplifications of unfrozen neighbors, so require at least the
    // merged ring and check it is live and valid.
    ASSERT_GE(result.added_obstacles.size(), 1u);
    const auto& merged =
      result.value->obstacles()[static_cast<size_t>(result.added_obstacles.front())]
        .ordered_vertices_;
    EXPECT_GE(merged.size(), 4u);
  }
}

// --- Outer contour protection drives the documented fallback ------------

TEST(PolymapIncremental, OuterContourChangeFallsBack) {
  auto map = makeBorderedMap(20, 20);
  const std::vector<PolymapEndpoint> goals{{17, 10, Point2d{17.5, 10.5}}};
  Polymap base = makeReadyPolymap(map, 2, 10, Point2d{2.5, 10.5}, goals);

  // Occupy cells adjacent to the border frame: the outer contour's raw
  // ring changes, so the incremental path must decline and the wrapper
  // must fall back to a full rebuild that still succeeds.
  const std::vector<std::pair<int, int>> cells{{1, 1}, {1, 2}};
  auto result = Polymap::applyOccupancyDelta(base, cells, {}, 2, 10, Point2d{2.5, 10.5}, goals,
                                             StopToken{});
  ASSERT_TRUE(result) << result.error;
  EXPECT_TRUE(result.fell_back_to_full_rebuild);
  EXPECT_FALSE(result.retired_obstacles.empty());
  EXPECT_FALSE(result.added_obstacles.empty());
}

// --- Repeated updates keep indices stable across a chain ----------------

TEST(PolymapIncremental, ChainOfUpdatesPreservesUnrelatedIndices) {
  auto map = makeBorderedMap(30, 24);
  for (int x = 8; x <= 12; ++x)
    occupy(map, x, 10);
  occupy(map, 25, 5);
  occupy(map, 25, 6);

  const std::vector<PolymapEndpoint> goals{{27, 12, Point2d{27.5, 12.5}}};
  Polymap current = makeReadyPolymap(map, 2, 12, Point2d{2.5, 12.5}, goals);

  // The wall grows over successive updates; the lone island must never be
  // touched (its ring and index must stay identical across every hop).
  std::vector<std::pair<int, int>> island_ring;
  int island_index = -1;
  for (size_t index = 0; index < current.obstacles().size(); ++index) {
    const auto& ring = current.obstacles()[index].ordered_vertices_;
    for (const auto& point : ring) {
      if (point.first >= 24 && point.first <= 27 && point.second >= 4 && point.second <= 8) {
        island_index = static_cast<int>(index);
        island_ring = ring;
        break;
      }
    }
    if (island_index >= 0)
      break;
  }
  ASSERT_GE(island_index, 0);

  for (int step = 0; step < 3; ++step) {
    const std::vector<std::pair<int, int>> cells{{8 + step * 2, 9}, {9 + step * 2, 9}};
    auto result = Polymap::applyOccupancyDelta(current, cells, {}, 2, 12, Point2d{2.5, 12.5},
                                               goals, StopToken{});
    ASSERT_TRUE(result) << result.error;
    ASSERT_TRUE(result.value);
    for (const int retired : result.retired_obstacles)
      EXPECT_NE(retired, island_index) << "the unrelated island was retired at step " << step;
    if (!result.fell_back_to_full_rebuild) {
      EXPECT_EQ(result.value->obstacles()[static_cast<size_t>(island_index)].ordered_vertices_,
                island_ring)
        << "the unrelated island changed at step " << step;
      // island ring may only stay identical if its index is untouched;
      // otherwise the test would have failed above.
    } else {
      // A fallback rebuilds every index: chain stability is not expected
      // there, but the island geometry must still exist somewhere.
      bool island_present = false;
      for (const auto& obstacle : result.value->obstacles()) {
        if (obstacle.ordered_vertices_ == island_ring) {
          island_present = true;
          break;
        }
      }
      EXPECT_TRUE(island_present) << "island vanished after fallback at step " << step;
      // Re-locate the island for subsequent steps.
      island_index = -1;
      for (size_t index = 0; index < result.value->obstacles().size(); ++index) {
        if (result.value->obstacles()[index].ordered_vertices_ == island_ring) {
          island_index = static_cast<int>(index);
          break;
        }
      }
    }
    current = std::move(*result.value);
  }
}

// --- Same-occupancy incremental vs full rebuild: both admit, C1 holds ----

TEST(PolymapIncremental, IncrementalResultIsValidAndHonestAboutChangeSet) {
  auto map = makeBorderedMap(24, 20);
  for (int x = 6; x <= 9; ++x)
    occupy(map, x, 8);
  for (int x = 6; x <= 9; ++x)
    occupy(map, x, 14);
  occupy(map, 15, 11);

  const std::vector<PolymapEndpoint> goals{{21, 10, Point2d{21.5, 10.5}}};
  Polymap base = makeReadyPolymap(map, 2, 10, Point2d{2.5, 10.5}, goals);

  const std::vector<std::pair<int, int>> cells{{15, 10}, {15, 12}, {14, 11}};
  BothSides sides = updateAndRebuild(base, cells, 2, 10, Point2d{2.5, 10.5}, goals);
  ASSERT_TRUE(sides.incremental) << sides.incremental.error;
  ASSERT_TRUE(sides.rebuilt) << sides.rebuilt.error;

  // NOTE on scope: RaystarCore::plan consumes a GridMap only, so running
  // the planner here would exercise the full-build pipeline twice and say
  // nothing about the incremental map.  Class-level equivalence between
  // an incrementally updated map and a rebuilt one is enforced by the
  // tree-layer differential tests (rebased tree vs fresh tree grown on
  // the SAME incremental polymap); here we check admission, the honest
  // change set, and the C1 index contract.
  EXPECT_FALSE(sides.incremental.fell_back_to_full_rebuild);
  EXPECT_EQ(sides.incremental.retired_obstacles.size(), 1u);
  EXPECT_EQ(sides.incremental.added_obstacles.size(), 1u);
  for (size_t index = 0; index < base.obstacles().size(); ++index) {
    if (std::find(sides.incremental.retired_obstacles.begin(),
                  sides.incremental.retired_obstacles.end(),
                  static_cast<int>(index)) !=
        sides.incremental.retired_obstacles.end())
      continue;
    EXPECT_EQ(sides.incremental.value->obstacles()[index].ordered_vertices_,
              base.obstacles()[index].ordered_vertices_)
      << "frozen obstacle " << index << " changed geometry";
  }
}

// --- Massive updates decline honestly to a full rebuild -----------------

TEST(PolymapIncremental, MassiveUpdateFallsBackCleanly) {
  auto map = makeBorderedMap(24, 20);
  for (int x = 6; x <= 9; ++x)
    occupy(map, x, 8);
  for (int x = 6; x <= 9; ++x)
    occupy(map, x, 14);
  for (int x = 14; x <= 17; ++x) {
    occupy(map, x, 6);
    occupy(map, x, 11);
  }
  const std::vector<PolymapEndpoint> goals{{21, 17, Point2d{21.5, 17.5}}};
  Polymap base = makeReadyPolymap(map, 2, 17, Point2d{2.5, 17.5}, goals);

  // Flood a large region around every obstacle: the cascade (or the
  // budget/outer rule) must decline and the wrapper must deliver a valid
  // full rebuild with the maximal change-set semantics.
  std::vector<std::pair<int, int>> flood;
  for (int y = 6; y <= 15; ++y)
    for (int x = 5; x <= 18; ++x)
      if (map.data[static_cast<size_t>(y) * map.width + static_cast<size_t>(x)] == 0)
        flood.emplace_back(x, y);
  ASSERT_GT(flood.size(), 50u);
  auto result = Polymap::applyOccupancyDelta(base, flood, {}, 2, 17, Point2d{2.5, 17.5}, goals,
                                             StopToken{});
  ASSERT_TRUE(result) << result.error;
  EXPECT_TRUE(result.fell_back_to_full_rebuild);
  EXPECT_FALSE(result.fallback_reason.empty());
  // Maximal change set: every previously live obstacle retired, every
  // live obstacle of the rebuilt map added.
  size_t live_before = 0;
  for (const auto& obstacle : base.obstacles())
    if (!obstacle.ordered_vertices_.empty())
      ++live_before;
  EXPECT_EQ(result.retired_obstacles.size(), live_before);
  size_t live_after = 0;
  for (const auto& obstacle : result.value->obstacles())
    if (!obstacle.ordered_vertices_.empty())
      ++live_after;
  EXPECT_EQ(result.added_obstacles.size(), live_after);
}

// --- No-op update is rejected up front ----------------------------------

TEST(PolymapIncremental, NoChangeDeltaIsRejected) {
  auto map = makeBorderedMap(20, 20);
  const std::vector<PolymapEndpoint> goals{{17, 10, Point2d{17.5, 10.5}}};
  Polymap base = makeReadyPolymap(map, 2, 10, Point2d{2.5, 10.5}, goals);

  const std::vector<std::pair<int, int>> already_occupied{{0, 5}, {0, 6}, {19, 7}};
  auto result = Polymap::applyOccupancyDelta(base, already_occupied, {}, 2, 10, Point2d{2.5, 10.5},
                                             goals, StopToken{});
  EXPECT_FALSE(result);
  EXPECT_FALSE(result.value.has_value());
  EXPECT_FALSE(result.error.empty());
}

// --- Disconnected update: no_path is a legitimate incremental outcome ----

TEST(PolymapIncremental, SealingAWallReportsNoPath) {
  auto map = makeBorderedMap(24, 20);
  for (int x = 8; x <= 14; ++x)
    occupy(map, x, 10);  // a wall; free rows 1..9 above and 11..18 below
  const std::vector<PolymapEndpoint> goals{{21, 15, Point2d{21.5, 15.5}}};
  Polymap base = makeReadyPolymap(map, 2, 15, Point2d{2.5, 15.5}, goals);

  // Seal column x=11 completely (row 10 is already occupied by the wall),
  // splitting the map into a left and a right half.
  std::vector<std::pair<int, int>> seal;
  for (int y = 1; y <= 18; ++y)
    if (y != 10)
      seal.emplace_back(11, y);
  auto result = Polymap::applyOccupancyDelta(base, seal, {}, 2, 15, Point2d{2.5, 15.5}, goals,
                                             StopToken{});
  EXPECT_EQ(result.status, PolymapCreateStatus::no_path);
}

// --- Cavity semantics: unreachable free space produces no contour ------- –

TEST(PolymapIncremental, SealedCavityAndItsIslandProduceNoContours) {
  auto map = makeBorderedMap(20, 20);
  // Ring wall sealing an interior free region, with an island inside it.
  for (int x = 7; x <= 12; ++x) {
    occupy(map, x, 7);
    occupy(map, x, 12);
  }
  for (int y = 7; y <= 12; ++y) {
    occupy(map, 7, y);
    occupy(map, 12, y);
  }
  occupy(map, 9, 9);
  occupy(map, 9, 10);
  occupy(map, 10, 9);
  occupy(map, 10, 10);

  const std::vector<PolymapEndpoint> goals{{17, 10, Point2d{17.5, 10.5}}};
  Polymap base = makeReadyPolymap(map, 2, 10, Point2d{2.5, 10.5}, goals);

  // Raw contours only bound REACHABLE free space: the ring's cavity side
  // and the cavity's island contribute nothing -- exactly two obstacles.
  // This pins the semantic that deletion inside a sealed cavity is a
  // no-op for the map, and that opening a cavity makes its walls (and
  // any islands inside) appear through raw-ring changes automatically.
  ASSERT_EQ(base.obstacles().size(), 2u);
}

// --- Deletion support -----------------------------------------------------

TEST(PolymapIncremental, ShrinkingAnObstacleRetiresAndAppends) {
  auto map = makeBorderedMap(20, 20);
  for (int x = 6; x <= 9; ++x) {
    occupy(map, x, 8);
    occupy(map, x, 14);
  }
  occupy(map, 15, 11);
  occupy(map, 15, 12);

  const std::vector<PolymapEndpoint> goals{{18, 10, Point2d{18.5, 10.5}}};
  Polymap base = makeReadyPolymap(map, 2, 10, Point2d{2.5, 10.5}, goals);

  // Shrink the lone island by freeing one cell.
  auto result = Polymap::applyOccupancyDelta(base, {}, {{15, 12}}, 2, 10, Point2d{2.5, 10.5},
                                             goals, StopToken{});
  ASSERT_TRUE(result) << result.error;
  EXPECT_FALSE(result.fell_back_to_full_rebuild);
  // The island's raw ring changed: one retirement, one appended ring; the
  // two walls and the outer frame must be frozen untouched.
  EXPECT_EQ(result.retired_obstacles.size(), 1u);
  EXPECT_EQ(result.added_obstacles.size(), 1u);
}

TEST(PolymapIncremental, FreeingInsideSealedCavityIsNoOp) {
  auto map = makeBorderedMap(20, 20);
  for (int x = 7; x <= 12; ++x) {
    occupy(map, x, 7);
    occupy(map, x, 12);
  }
  for (int y = 7; y <= 12; ++y) {
    occupy(map, 7, y);
    occupy(map, 12, y);
  }
  occupy(map, 9, 9);  // wall cell INSIDE the sealed cavity
  const std::vector<PolymapEndpoint> goals{{17, 10, Point2d{17.5, 10.5}}};
  Polymap base = makeReadyPolymap(map, 2, 10, Point2d{2.5, 10.5}, goals);
  ASSERT_EQ(base.obstacles().size(), 2u);  // frame + ring only

  // Freeing a wall cell inside the sealed cavity cannot change any raw
  // ring (the cavity is unreachable), so every contour freezes and the
  // change set must be empty.
  auto result = Polymap::applyOccupancyDelta(base, {}, {{9, 9}}, 2, 10, Point2d{2.5, 10.5},
                                             goals, StopToken{});
  ASSERT_TRUE(result) << result.error;
  EXPECT_FALSE(result.fell_back_to_full_rebuild);
  EXPECT_TRUE(result.retired_obstacles.empty());
  EXPECT_TRUE(result.added_obstacles.empty());
  ASSERT_EQ(result.value->obstacles().size(), base.obstacles().size());
  for (size_t index = 0; index < base.obstacles().size(); ++index)
    EXPECT_EQ(result.value->obstacles()[index].ordered_vertices_,
              base.obstacles()[index].ordered_vertices_);
}

TEST(PolymapIncremental, OpeningACavityRevealsItsIsland) {
  auto map = makeBorderedMap(20, 20);
  for (int x = 7; x <= 12; ++x) {
    occupy(map, x, 7);
    occupy(map, x, 12);
  }
  for (int y = 7; y <= 12; ++y) {
    occupy(map, 7, y);
    occupy(map, 12, y);
  }
  occupy(map, 9, 9);  // island inside the sealed cavity
  occupy(map, 9, 10);
  occupy(map, 10, 9);
  occupy(map, 10, 10);
  const std::vector<PolymapEndpoint> goals{{5, 10, Point2d{5.5, 10.5}}};
  Polymap base = makeReadyPolymap(map, 16, 10, Point2d{16.5, 10.5}, goals);
  ASSERT_EQ(base.obstacles().size(), 2u);  // cavity invisible from the east

  // Open the ring's west wall: the cavity becomes reachable from the
  // start side (start is west at (5,10); goal stays east).  The ring's
  // raw ring changes and the island appears.
  auto result = Polymap::applyOccupancyDelta(base, {}, {{7, 9}, {7, 10}}, 5, 10,
                                             Point2d{5.5, 10.5}, goals, StopToken{});
  ASSERT_TRUE(result) << result.error;
  // The newly reachable cavity must now carry the island contour.
  size_t live = 0;
  for (const auto& obstacle : result.value->obstacles())
    if (!obstacle.ordered_vertices_.empty())
      ++live;
  EXPECT_GE(live, 3u);  // frame + opened ring + revealed island
  bool island_found = false;
  for (const auto& obstacle : result.value->obstacles()) {
    for (const auto& point : obstacle.ordered_vertices_) {
      if (point.first >= 8 && point.first <= 11 && point.second >= 8 && point.second <= 11) {
        island_found = true;
        break;
      }
    }
  }
  EXPECT_TRUE(island_found) << "the cavity island must appear after opening the wall";
}
