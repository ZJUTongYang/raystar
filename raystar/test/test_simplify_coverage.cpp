// Rasterized conservative-cover property test (SIMPLIFY_STACK_SCAN_DESIGN.md
// section 6.2): after Polymap construction -- extraction plus greedy
// simplification -- every occupied cell of the working occupancy must still
// be covered by the polygon set, and the protected start/goal positions must
// lie strictly in free space.  This is the direct per-instance check of
// R-L1's conservative cover: it depends only on the OUTPUT, not on the
// removal order, so it guards any future simplification change.
//
// Coverage semantics: rings are free-space boundaries.  The outer ring's
// bounded side is the free interior, interior rings' bounded sides are
// obstacles.  An occupied cell point is covered iff it is inside-or-on some
// interior ring, or not strictly inside the outer ring.  Cell centres and
// all four corners are sampled; boundaries count as obstacle closure.

#include <gtest/gtest.h>
#include <raystar/polymap.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>

#include <cstdint>
#include <utility>
#include <vector>

#include "simplify_test_worlds.h"

namespace raystar {

// Friend peer (polymap.h): re-runs the private simplification entry point
// for the maximality assertion.
class PolymapTestPeer {
public:
  static void resimplify(Polymap& polymap, const Point2d& start, const Point2d& goal) {
    polymap.simplifyPolyObstacles(start, goal);
  }
};

namespace {

using simplify_test_worlds::connected;
using simplify_test_worlds::randomWorld;
using simplify_test_worlds::Rng;

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Polygon = CGAL::Polygon_2<Kernel>;

TEST(SimplifyCoverage, EveryOccupiedCellStaysCoveredAndEndpointsStayFree) {
  int worlds_checked = 0;
  for (std::uint64_t seed = 1; seed <= 40; ++seed) {
    Rng rng{20260926ULL * 31 + seed};
    const int width = rng.range(36, 60);
    const int height = rng.range(30, 48);
    const GridMap map = randomWorld(rng, width, height);

    // A connected free start/goal pair (cell centres).
    int start_x = -1, start_y = -1, goal_x = -1, goal_y = -1;
    for (int attempt = 0; attempt < 60; ++attempt) {
      const int sx = rng.range(1, width - 2), sy = rng.range(1, height - 2);
      const int gx = rng.range(1, width - 2), gy = rng.range(1, height - 2);
      if ((sx == gx && sy == gy) ||
          map.data[static_cast<size_t>(sy) * map.width + static_cast<size_t>(sx)] != 0 ||
          map.data[static_cast<size_t>(gy) * map.width + static_cast<size_t>(gx)] != 0)
        continue;
      if (connected(map, sx, sy, gx, gy)) {
        start_x = sx;
        start_y = sy;
        goal_x = gx;
        goal_y = gy;
        break;
      }
    }
    if (start_x < 0)
      continue;
    const Point2d start{start_x + 0.5, start_y + 0.5};
    const Point2d goal{goal_x + 0.5, goal_y + 0.5};

    auto created = Polymap::create(
      map, start_x, start_y, start, {PolymapEndpoint{goal_x, goal_y, goal}}, StopToken{});
    if (!created)
      continue;  // degenerate world rejected upstream: nothing to assert
    const Polymap& polymap = *created.value;
    ++worlds_checked;

    // The construction outlines the border; assert against ITS occupancy.
    const auto& occupancy = polymap.occupancyData();
    const int poly_width = polymap.width();
    const int poly_height = polymap.height();

    // Rings as exact polygons; the outer ring is the one whose bounded side
    // contains the start (free space lives inside it).
    std::vector<Polygon> rings;
    int outer = -1;
    for (const auto& obstacle : polymap.obstacles()) {
      Polygon polygon;
      for (const auto& vertex : obstacle.ordered_vertices_)
        polygon.push_back(Kernel::Point_2(vertex.first, vertex.second));
      if (polygon.size() >= 3 &&
          polygon.bounded_side(Kernel::Point_2(start.first, start.second)) == CGAL::ON_BOUNDED_SIDE)
        outer = static_cast<int>(rings.size());
      rings.push_back(std::move(polygon));
    }
    ASSERT_GE(outer, 0) << "seed " << seed << ": no outer ring contains the start";

    const auto covered = [&](double x, double y) {
      const Kernel::Point_2 point(x, y);
      for (size_t ring = 0; ring < rings.size(); ++ring) {
        if (rings[ring].size() < 3)
          continue;
        const auto side = rings[ring].bounded_side(point);
        if (static_cast<int>(ring) == outer) {
          if (side != CGAL::ON_BOUNDED_SIDE)
            return true;  // on or outside the free-space boundary
        } else if (side != CGAL::ON_UNBOUNDED_SIDE) {
          return true;  // inside or on an interior obstacle polygon
        }
      }
      return false;
    };

    for (int y = 0; y < poly_height; ++y) {
      for (int x = 0; x < poly_width; ++x) {
        if (occupancy[static_cast<size_t>(y) * static_cast<size_t>(poly_width) +
                      static_cast<size_t>(x)] == 0)
          continue;
        EXPECT_TRUE(covered(x + 0.5, y + 0.5))
          << "seed " << seed << ": centre of occupied cell (" << x << "," << y
          << ") escaped the simplified polygons";
        for (const auto& corner :
             {std::pair<int, int>{x, y}, {x + 1, y}, {x, y + 1}, {x + 1, y + 1}}) {
          EXPECT_TRUE(covered(corner.first, corner.second))
            << "seed " << seed << ": corner (" << corner.first << "," << corner.second
            << ") of occupied cell (" << x << "," << y << ") escaped";
        }
      }
    }

    // Protected endpoints must remain strictly in free space: strictly
    // inside the outer ring and strictly outside every interior ring.
    for (const auto& endpoint : {start, goal}) {
      const Kernel::Point_2 point(endpoint.first, endpoint.second);
      EXPECT_EQ(rings[static_cast<size_t>(outer)].bounded_side(point), CGAL::ON_BOUNDED_SIDE)
        << "seed " << seed;
      for (size_t ring = 0; ring < rings.size(); ++ring) {
        if (static_cast<int>(ring) == outer || rings[ring].size() < 3)
          continue;
        EXPECT_EQ(rings[ring].bounded_side(point), CGAL::ON_UNBOUNDED_SIDE)
          << "seed " << seed << " ring " << ring;
      }
    }
  }
  // The property must have been exercised on a healthy share of the seeds.
  EXPECT_GE(worlds_checked, 20);
}

TEST(SimplifyCoverage, SecondSimplifyPassRemovesNothing) {
  // Greedy maximality (design Proposition 1): the loop terminates only when
  // a full lap finds nothing removable, so re-simplifying the output with
  // the same protected points must be a no-op.  A broken confirmation lap
  // (early exit, missed far unlock) fails this on some seed.
  for (std::uint64_t seed = 1; seed <= 15; ++seed) {
    Rng rng{20260926ULL * 13 + seed};
    const int width = rng.range(36, 56);
    const int height = rng.range(30, 44);
    const GridMap map = randomWorld(rng, width, height);
    int start_x = -1, start_y = -1, goal_x = -1, goal_y = -1;
    for (int attempt = 0; attempt < 60 && start_x < 0; ++attempt) {
      const int sx = rng.range(1, width - 2), sy = rng.range(1, height - 2);
      const int gx = rng.range(1, width - 2), gy = rng.range(1, height - 2);
      if ((sx == gx && sy == gy) ||
          map.data[static_cast<size_t>(sy) * map.width + static_cast<size_t>(sx)] != 0 ||
          map.data[static_cast<size_t>(gy) * map.width + static_cast<size_t>(gx)] != 0)
        continue;
      if (connected(map, sx, sy, gx, gy)) {
        start_x = sx;
        start_y = sy;
        goal_x = gx;
        goal_y = gy;
      }
    }
    if (start_x < 0)
      continue;
    const Point2d start{start_x + 0.5, start_y + 0.5};
    const Point2d goal{goal_x + 0.5, goal_y + 0.5};
    auto created = Polymap::create(
      map, start_x, start_y, start, {PolymapEndpoint{goal_x, goal_y, goal}}, StopToken{});
    if (!created)
      continue;
    Polymap& polymap = *created.value;
    size_t before = 0;
    for (const auto& obstacle : polymap.obstacles()) before += obstacle.ordered_vertices_.size();
    PolymapTestPeer::resimplify(polymap, start, goal);
    size_t after = 0;
    for (const auto& obstacle : polymap.obstacles()) after += obstacle.ordered_vertices_.size();
    EXPECT_EQ(before, after) << "seed " << seed << ": the first pass was not greedy-maximal";
  }
}

}  // namespace
}  // namespace raystar
