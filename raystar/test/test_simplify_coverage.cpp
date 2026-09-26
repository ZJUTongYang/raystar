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

namespace raystar {
namespace {

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Polygon = CGAL::Polygon_2<Kernel>;

struct Rng {
  std::uint64_t state;
  std::uint32_t next() {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<std::uint32_t>(state >> 33);
  }
  int range(int lo, int hi) {  // inclusive
    return lo + static_cast<int>(next() % static_cast<std::uint32_t>(hi - lo + 1));
  }
};

GridMap randomWorld(Rng& rng, int width, int height) {
  GridMap map;
  map.width = static_cast<unsigned int>(width);
  map.height = static_cast<unsigned int>(height);
  map.resolution = 1.0F;
  map.data.assign(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
  const auto block = [&](int x0, int y0, int x1, int y1) {
    for (int y = y0; y <= y1; ++y)
      for (int x = x0; x <= x1; ++x)
        if (x >= 0 && y >= 0 && x < width && y < height)
          map.data[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)] =
            1;
  };
  block(0, 0, width - 1, 0);
  block(0, height - 1, width - 1, height - 1);
  block(0, 0, 0, height - 1);
  block(width - 1, 0, width - 1, height - 1);
  const int walls = rng.range(2, 5);
  for (int wall = 0; wall < walls; ++wall) {
    if (rng.range(0, 1) == 1) {
      const int y = rng.range(4, height - 5);
      const int x0 = rng.range(2, width - 12);
      block(x0, y, x0 + rng.range(6, 14), y);
    } else {
      const int x = rng.range(4, width - 5);
      const int y0 = rng.range(2, height - 12);
      block(x, y0, x, y0 + rng.range(6, 14));
    }
  }
  const int blobs = rng.range(1, 4);
  for (int blob = 0; blob < blobs; ++blob) {
    const int x = rng.range(2, width - 6);
    const int y = rng.range(2, height - 6);
    block(x, y, x + rng.range(1, 3), y + rng.range(1, 3));
  }
  return map;
}

bool connected(const GridMap& map, int sx, int sy, int gx, int gy) {
  const int width = static_cast<int>(map.width);
  const int height = static_cast<int>(map.height);
  std::vector<char> seen(map.data.size(), 0);
  std::vector<int> stack{sy * width + sx};
  seen[static_cast<size_t>(sy * width + sx)] = 1;
  while (!stack.empty()) {
    const int cell = stack.back();
    stack.pop_back();
    if (cell == gy * width + gx)
      return true;
    const int x = cell % width;
    const int y = cell / width;
    const int neighbours[4][2] = {{x - 1, y}, {x + 1, y}, {x, y - 1}, {x, y + 1}};
    for (const auto& neighbour : neighbours) {
      if (neighbour[0] < 0 || neighbour[1] < 0 || neighbour[0] >= width || neighbour[1] >= height)
        continue;
      const int index = neighbour[1] * width + neighbour[0];
      if (map.data[static_cast<size_t>(index)] != 0 || seen[static_cast<size_t>(index)])
        continue;
      seen[static_cast<size_t>(index)] = 1;
      stack.push_back(index);
    }
  }
  return false;
}

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

}  // namespace
}  // namespace raystar
