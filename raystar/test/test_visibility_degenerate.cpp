#include <gtest/gtest.h>
#include <raystar/exact_geometry.h>
#include <raystar/polymap.h>
#include <raystar/raystar_core.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace {

using raystar::exactPoint;
using raystar::GridMap;
using raystar::PlanResult;
using raystar::Point2d;
using raystar::Polymap;
using raystar::PolymapCreateResult;
using raystar::RaystarCore;
using raystar::VisibilityRegion;

// The same regularization cnlspp applies to its integer endpoints: the
// x:y ratio 1:sqrt(2) is irrational, so the shift breaks collinearity with
// every rational-slope line through integer grid vertices.  The perturbed
// twin is therefore a non-degenerate reference for what the exact integer
// source must see by continuity.
constexpr double kPerturbX = 1e-6;
constexpr double kPerturbY = 1e-6 * 1.4142135623730951;

GridMap makeBorderedMap(unsigned int width, unsigned int height) {
  GridMap map;
  map.width = width;
  map.height = height;
  map.resolution = 1.0f;
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

void occupyCell(GridMap& map, unsigned int x, unsigned int y) {
  map.data[static_cast<size_t>(y) * map.width + x] = 1;
}

// +x ray from the integer source s=(2,5) lies exactly on the top edges of
// two obstacles below it: a=(4,5)-b=(7,5) and c=(8,5)-d=(10,5), with the
// free gap b-c on the ray between them.
GridMap makeTangentGapMap() {
  auto map = makeBorderedMap(12, 12);
  for (unsigned int x = 4; x <= 6; ++x)
    occupyCell(map, x, 4);
  occupyCell(map, 8, 4);
  occupyCell(map, 9, 4);
  return map;
}

// Sandwich variant: the first grazed edge lies below the ray and the second
// above it, so the ray threads the pinch between the two obstacles.
GridMap makeTangentSandwichMap() {
  auto map = makeBorderedMap(12, 12);
  for (unsigned int x = 4; x <= 6; ++x)
    occupyCell(map, x, 4);
  occupyCell(map, 8, 5);
  occupyCell(map, 9, 5);
  return map;
}

// Corridor variant: the ray runs along the top edge of a long wall with a
// parallel wall above it.
GridMap makeCorridorFloorMap() {
  auto map = makeBorderedMap(12, 12);
  for (unsigned int x = 3; x <= 8; ++x) {
    occupyCell(map, x, 4);
    occupyCell(map, x, 6);
  }
  return map;
}

GridMap makeSingleTangentEdgeMap() {
  auto map = makeBorderedMap(12, 12);
  occupyCell(map, 5, 4);
  occupyCell(map, 6, 4);
  return map;
}

// Reproduction of the failing random batch map (seed 20260918, map index 6):
// the +x ray from the integer source s=(6,4) grazes the BOTTOM edge of the
// 2x2 block from the free side below -- the mirror of the configurations
// above, and the one that currently makes the sweep fail closed.
GridMap makeBottomGrazedBlockMap() {
  auto map = makeBorderedMap(12, 12);
  occupyCell(map, 7, 2);
  occupyCell(map, 8, 2);
  occupyCell(map, 7, 3);
  occupyCell(map, 8, 3);
  occupyCell(map, 5, 7);
  occupyCell(map, 5, 8);
  return map;
}

// The caller keeps the result (and thereby the non-copyable Polymap) alive;
// planning mutates the map's caches, so a non-const reference is required.
PolymapCreateResult createPolymapForRoot(const GridMap& map, const Point2d& source) {
  const auto start_cell_x = static_cast<int>(std::floor(source.first));
  const auto start_cell_y = static_cast<int>(std::floor(source.second));
  return Polymap::create(map, start_cell_x, start_cell_y, 2, 8, source, Point2d{2.5, 8.5});
}

int findEndpoint(const VisibilityRegion& region, double x, double y) {
  for (size_t i = 0; i < region.size(); ++i) {
    if (region[i].position.first == x && region[i].position.second == y)
      return static_cast<int>(i);
  }
  return -1;
}

bool regionContains(const VisibilityRegion& region, double x, double y) {
  return findEndpoint(region, x, y) >= 0;
}

void expectRootRegionValidates(Polymap& polymap, const Point2d& source) {
  VisibilityRegion region;
  std::string error;
  ASSERT_TRUE(polymap.getRootVisibilityRegion(source, region, &error)) << error;
  std::string validation_error;
  EXPECT_TRUE(polymap.validateVisibilityRegion(region, &validation_error)) << validation_error;
}

// The exact integer source must see the same homotopy-class representatives
// as its infinitesimally perturbed twin: identical feasibility, identical
// class count, and matching certified costs (geometry shifts by O(epsilon)).
void expectTopKMatchesPerturbedTwin(const GridMap& map,
                                    const Point2d& source,
                                    const Point2d& goal,
                                    int k) {
  const Point2d shifted{source.first + kPerturbX, source.second + kPerturbY};
  RaystarCore core;
  const PlanResult exact = core.plan(map, source, goal, k, false);
  const PlanResult shifted_result = core.plan(map, shifted, goal, k, false);

  ASSERT_EQ(exact.success, shifted_result.success)
    << "feasibility differs between the exact and perturbed source"
    << "; exact message: " << exact.message
    << "; shifted message: " << shifted_result.message;
  if (!exact.success)
    return;
  ASSERT_EQ(exact.path_solutions.size(), shifted_result.path_solutions.size())
    << "homotopy-class count differs between the exact and perturbed source";
  for (size_t i = 0; i < exact.path_solutions.size(); ++i) {
    EXPECT_NEAR(exact.path_solutions[i].path_cost_,
                shifted_result.path_solutions[i].path_cost_,
                // The twin's geometry shifts by O(|perturbation|) ~ 1.7e-6,
                // so allow several perturbation lengths; a missing or extra
                // homotopy class changes costs at the 1e-1..1e0 scale.
                1e-5 * (1.0 + exact.path_solutions[i].path_cost_))
      << "class " << i << " cost differs between the exact and perturbed source";
    // The reported start stays within the general-position recovery offset
    // of the requested source (equal to it when no recovery was needed).
    EXPECT_NEAR(exact.path_solutions[i].start_.first, source.first, 2e-6);
    EXPECT_NEAR(exact.path_solutions[i].start_.second, source.second, 2e-6);
  }
}

struct Lcg {
  std::uint64_t state;
  std::uint32_t next() {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<std::uint32_t>(state >> 33);
  }
};

// Random bordered map with axis-aligned integer rectangles plus an integer
// source strictly inside free space and a cell-centre goal.  Axis-aligned
// rectangles make tangent rays from the integer source common, which is the
// degenerate family this suite guards.
bool makeRandomTangentMap(Lcg& rng,
                          GridMap& map,
                          Point2d& source,
                          Point2d& goal) {
  map = makeBorderedMap(12, 12);
  const int rect_count = 3 + static_cast<int>(rng.next() % 3);
  for (int rect = 0; rect < rect_count; ++rect) {
    const unsigned int width = 1 + rng.next() % 2;
    const unsigned int height = 1 + rng.next() % 2;
    const unsigned int x0 = 2 + rng.next() % 8;
    const unsigned int y0 = 2 + rng.next() % 8;
    for (unsigned int y = y0; y < std::min(y0 + height, 11u); ++y)
      for (unsigned int x = x0; x < std::min(x0 + width, 11u); ++x)
        occupyCell(map, x, y);
  }

  const auto cells_free_around = [&map](unsigned int x, unsigned int y) {
    return map.data[(y - 1) * map.width + (x - 1)] == 0 &&
           map.data[(y - 1) * map.width + x] == 0 &&
           map.data[y * map.width + (x - 1)] == 0 && map.data[y * map.width + x] == 0;
  };

  bool have_source = false;
  for (int attempt = 0; attempt < 64 && !have_source; ++attempt) {
    const unsigned int x = 2 + rng.next() % 8;
    const unsigned int y = 2 + rng.next() % 8;
    if (cells_free_around(x, y)) {
      source = Point2d{static_cast<double>(x), static_cast<double>(y)};
      have_source = true;
    }
  }
  if (!have_source)
    return false;

  for (int attempt = 0; attempt < 64; ++attempt) {
    const unsigned int x = 2 + rng.next() % 8;
    const unsigned int y = 2 + rng.next() % 8;
    if (map.data[y * map.width + x] == 0) {
      goal = Point2d{x + 0.5, y + 0.5};
      return true;
    }
  }
  return false;
}

}  // namespace

// --- Crafted tangent configurations: root visibility must succeed -------- –

TEST(TangentRayVisibility, GapBetweenCollinearEdgesIsPreserved) {
  const auto map = makeTangentGapMap();
  auto polymap_result = createPolymapForRoot(map, Point2d{2, 5});
  ASSERT_TRUE(polymap_result) << polymap_result.error;
  Polymap& polymap = *polymap_result.value;

  VisibilityRegion region;
  std::string error;
  ASSERT_TRUE(polymap.getRootVisibilityRegion(Point2d{2, 5}, region, &error)) << error;

  std::string validation_error;
  EXPECT_TRUE(polymap.validateVisibilityRegion(region, &validation_error)) << validation_error;

  EXPECT_TRUE(regionContains(region, 4, 5)) << "near corner a of the first grazed edge";
  EXPECT_TRUE(regionContains(region, 7, 5)) << "far corner b of the first grazed edge";
  EXPECT_TRUE(regionContains(region, 8, 5)) << "near corner c of the second grazed edge";
  EXPECT_TRUE(regionContains(region, 10, 5)) << "far corner d of the second grazed edge";

  // b and c must stay adjacent in the ordered region and on the same ray
  // from the source: they delimit the free gap the tree must be able to
  // expand through.  Collapsing the stacked events to a single cutoff point
  // would silently drop that homotopy branch.
  const raystar::exact_geometry::Point exact_source(2, 5);
  bool gap_adjacent = false;
  for (size_t i = 0; i + 1 < region.size(); ++i) {
    const bool forward = region[i].position.first == 7 && region[i].position.second == 5 &&
                         region[i + 1].position.first == 8 && region[i + 1].position.second == 5;
    const bool backward = region[i].position.first == 8 && region[i].position.second == 5 &&
                          region[i + 1].position.first == 7 && region[i + 1].position.second == 5;
    if ((forward || backward) &&
        raystar::exact_geometry::isSameDirectedRay(exact_source, exactPoint(region[i]),
                                                   exactPoint(region[i + 1]))) {
      gap_adjacent = true;
      break;
    }
  }
  EXPECT_TRUE(gap_adjacent)
    << "the on-ray gap b-c must remain a same-ray adjacent endpoint pair";
}

TEST(TangentRayVisibility, SandwichPinchSucceeds) {
  const auto map = makeTangentSandwichMap();
  auto polymap_result = createPolymapForRoot(map, Point2d{2, 5});
  ASSERT_TRUE(polymap_result) << polymap_result.error;
  expectRootRegionValidates(*polymap_result.value, Point2d{2, 5});
}

TEST(TangentRayVisibility, CorridorFloorEdgeIsFollowed) {
  const auto map = makeCorridorFloorMap();
  auto polymap_result = createPolymapForRoot(map, Point2d{2, 5});
  ASSERT_TRUE(polymap_result) << polymap_result.error;
  Polymap& polymap = *polymap_result.value;

  VisibilityRegion region;
  std::string error;
  ASSERT_TRUE(polymap.getRootVisibilityRegion(Point2d{2, 5}, region, &error)) << error;
  std::string validation_error;
  EXPECT_TRUE(polymap.validateVisibilityRegion(region, &validation_error)) << validation_error;
  // Both ends of the grazed floor edge must be part of the region.
  EXPECT_TRUE(regionContains(region, 3, 5));
  EXPECT_TRUE(regionContains(region, 9, 5));
}

TEST(TangentRayVisibility, SingleTangentEdgeSucceeds) {
  const auto map = makeSingleTangentEdgeMap();
  auto polymap_result = createPolymapForRoot(map, Point2d{2, 5});
  ASSERT_TRUE(polymap_result) << polymap_result.error;
  expectRootRegionValidates(*polymap_result.value, Point2d{2, 5});
}

// The direct Polymap-layer sweep deliberately stays fail-closed on this
// configuration (the exact same-ray group it produces has no consistent
// radial ordering).  RaystarCore::plan recovers through the
// general-position root retry, which the planner-level tests below verify.
// If the sweep itself is ever fixed to handle the tangency exactly, this
// marker test flips and should be updated together with the sweep fix.
TEST(TangentRayVisibility, BottomGrazedBlockDirectSweepFailsClosed) {
  const auto map = makeBottomGrazedBlockMap();
  // Goal matters: it is a simplification protection point, and the exact
  // tangency structure depends on which obstacle corners survive.  Use the
  // same goal as the planner-level reproduction below.
  auto polymap_result = Polymap::create(map, 6, 4, 6, 3, Point2d{6, 4}, Point2d{6.5, 3.5});
  ASSERT_TRUE(polymap_result) << polymap_result.error;
  Polymap& polymap = *polymap_result.value;

  VisibilityRegion region;
  std::string error;
  EXPECT_FALSE(polymap.getRootVisibilityRegion(Point2d{6, 4}, region, &error));
  EXPECT_TRUE(region.empty()) << "a failed sweep must clear the region";
}

// --- Planner-level differential: exact source vs perturbed twin -----------

TEST(TangentRayPlanning, TangentGapMapMatchesPerturbedTwin) {
  expectTopKMatchesPerturbedTwin(makeTangentGapMap(), Point2d{2, 5}, Point2d{2.5, 8.5}, 3);
}

TEST(TangentRayPlanning, SandwichMapMatchesPerturbedTwin) {
  expectTopKMatchesPerturbedTwin(makeTangentSandwichMap(), Point2d{2, 5}, Point2d{2.5, 8.5}, 3);
}

TEST(TangentRayPlanning, CorridorFloorMapMatchesPerturbedTwin) {
  expectTopKMatchesPerturbedTwin(makeCorridorFloorMap(), Point2d{2, 5}, Point2d{2.5, 8.5}, 3);
}

TEST(TangentRayPlanning, SingleTangentEdgeMapMatchesPerturbedTwin) {
  expectTopKMatchesPerturbedTwin(makeSingleTangentEdgeMap(), Point2d{2, 5}, Point2d{2.5, 8.5}, 3);
}

TEST(TangentRayPlanning, BottomGrazedBlockMapMatchesPerturbedTwin) {
  expectTopKMatchesPerturbedTwin(makeBottomGrazedBlockMap(), Point2d{6, 4}, Point2d{6.5, 3.5}, 3);
}

TEST(TangentRayPlanning, RandomAxisAlignedMapsMatchPerturbedTwin) {
  Lcg rng{20260918ULL};
  int checked = 0;
  for (int map_index = 0; map_index < 40; ++map_index) {
    GridMap map;
    Point2d source;
    Point2d goal;
    if (!makeRandomTangentMap(rng, map, source, goal))
      continue;
    ASSERT_NO_FATAL_FAILURE(expectTopKMatchesPerturbedTwin(map, source, goal, 3))
      << "random map " << map_index << " source (" << source.first << "," << source.second
      << ") goal (" << goal.first << "," << goal.second << ")";
    ++checked;
  }
  EXPECT_GE(checked, 20) << "the random batch should exercise a usable number of maps";
}
