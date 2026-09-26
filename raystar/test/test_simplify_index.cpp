// Unit tests for SimplifyCandidateIndex (SIMPLIFY_STACK_SCAN_DESIGN.md M1).
//
// The completeness lemma reduces the index's obligation to a specification:
// a query must return exactly the registrations whose bounding box overlaps
// the query box.  The differential tests below therefore compare the index
// against a brute-force bbox-overlap oracle on randomized ring sets; the
// remaining tests pin the removal bookkeeping and the indexability gate.

#include <gtest/gtest.h>
#include <raystar/polymap.h>

#include <algorithm>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "polymap_simplify_index.h"

namespace raystar {
namespace {

using polymap_impl::SimplifyCandidateIndex;

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

Obs makeRing(std::vector<std::pair<int, int>> vertices) {
  Obs obs;
  obs.ordered_vertices_ = std::move(vertices);
  return obs;
}

// Brute-force oracle: every vertex inside the closed box; every edge whose
// bbox overlaps the closed box.
void bruteForce(const std::vector<Obs>& obstacles,
                int min_x,
                int min_y,
                int max_x,
                int max_y,
                std::set<std::pair<int, int>>& vertices,
                std::set<std::pair<int, int>>& edges) {
  for (size_t obstacle = 0; obstacle < obstacles.size(); ++obstacle) {
    const auto& ring = obstacles[obstacle].ordered_vertices_;
    for (size_t vertex = 0; vertex < ring.size(); ++vertex) {
      const auto& coordinate = ring[vertex];
      if (coordinate.first >= min_x && coordinate.first <= max_x && coordinate.second >= min_y &&
          coordinate.second <= max_y)
        vertices.insert({static_cast<int>(obstacle), static_cast<int>(vertex)});
      const auto& to = ring[(vertex + 1) % ring.size()];
      const int edge_min_x = std::min(coordinate.first, to.first);
      const int edge_max_x = std::max(coordinate.first, to.first);
      const int edge_min_y = std::min(coordinate.second, to.second);
      const int edge_max_y = std::max(coordinate.second, to.second);
      if (edge_min_x <= max_x && min_x <= edge_max_x && edge_min_y <= max_y && min_y <= edge_max_y)
        edges.insert({static_cast<int>(obstacle), static_cast<int>(vertex)});
    }
  }
}

// Generates a random set of disjoint-coordinate rectangular-ish rings.
std::vector<Obs> randomRings(Rng& rng, int count) {
  std::vector<Obs> obstacles;
  std::set<std::pair<int, int>> used;
  for (int index = 0; index < count; ++index) {
    const int x = rng.range(0, 220);
    const int y = rng.range(0, 220);
    const int width = rng.range(2, 24);
    const int height = rng.range(2, 24);
    std::vector<std::pair<int, int>> ring = {
      {x, y}, {x + width, y}, {x + width, y + height}, {x, y + height}};
    // Optionally jitter one corner into a concavity (keeps coordinates
    // unique within the attempt).
    if (rng.range(0, 1) == 1)
      ring.push_back({x + width / 2, y + height / 2 + 1});
    bool clash = false;
    for (const auto& vertex : ring)
      if (used.count(vertex))
        clash = true;
    if (clash) {
      --index;
      continue;
    }
    for (const auto& vertex : ring) used.insert(vertex);
    obstacles.push_back(makeRing(std::move(ring)));
  }
  return obstacles;
}

TEST(SimplifyCandidateIndex, MatchesBruteForceOnRandomizedRingSets) {
  for (std::uint64_t seed = 1; seed <= 60; ++seed) {
    Rng rng{20260926ULL + seed};
    const auto obstacles = randomRings(rng, rng.range(2, 8));
    SimplifyCandidateIndex index(rng.range(1, 16));
    ASSERT_TRUE(index.build(obstacles)) << "seed " << seed;

    for (int query = 0; query < 40; ++query) {
      const int ax = rng.range(-10, 260), ay = rng.range(-10, 260);
      const int bx = rng.range(-10, 260), by = rng.range(-10, 260);
      const int min_x = std::min(ax, bx), max_x = std::max(ax, bx);
      const int min_y = std::min(ay, by), max_y = std::max(ay, by);

      std::set<std::pair<int, int>> expected_vertices, expected_edges;
      bruteForce(obstacles, min_x, min_y, max_x, max_y, expected_vertices, expected_edges);

      std::vector<SimplifyCandidateIndex::VertexRecord> vertices;
      std::vector<SimplifyCandidateIndex::EdgeRecord> edges;
      index.verticesInBox(min_x, min_y, max_x, max_y, vertices);
      index.edgesInBox(min_x, min_y, max_x, max_y, edges);

      std::set<std::pair<int, int>> got_vertices, got_edges;
      for (const auto& record : vertices) {
        EXPECT_TRUE(got_vertices.insert({record.obstacle, record.vertex}).second)
          << "duplicate vertex report, seed " << seed;
      }
      for (const auto& record : edges) {
        EXPECT_TRUE(got_edges.insert({record.obstacle, record.edge}).second)
          << "duplicate edge report, seed " << seed;
      }
      // Superset in the lemma direction is the correctness obligation; the
      // bucket discretization may only over-approximate.
      for (const auto& expected : expected_vertices)
        EXPECT_TRUE(got_vertices.count(expected))
          << "missing vertex (" << expected.first << "," << expected.second << ") seed " << seed;
      for (const auto& expected : expected_edges)
        EXPECT_TRUE(got_edges.count(expected))
          << "missing edge (" << expected.first << "," << expected.second << ") seed " << seed;
    }
  }
}

TEST(SimplifyCandidateIndex, DegenerateAndAdversarialQueries) {
  // A long diagonal edge, a vertical edge collinear with a query chord, and
  // an outer-frame ring corner.
  const std::vector<Obs> obstacles = {
    makeRing({{0, 0}, {100, 0}, {100, 100}, {0, 100}}),  // frame
    makeRing({{10, 20}, {90, 80}, {10, 80}}),            // long diagonal edge 0
    makeRing({{50, 5}, {52, 5}, {52, 15}, {50, 15}}),    // tall thin box
  };
  SimplifyCandidateIndex index(8);
  ASSERT_TRUE(index.build(obstacles));

  // Zero-area query (degenerate triangle => bbox of a chord): a horizontal
  // chord crossing the diagonal's midsection must report the diagonal edge.
  std::vector<SimplifyCandidateIndex::EdgeRecord> edges;
  index.edgesInBox(40, 50, 60, 50, edges);
  EXPECT_TRUE(std::any_of(edges.begin(), edges.end(), [](const auto& record) {
    return record.obstacle == 1 && record.edge == 0;
  }));

  // A chord collinear with and overlapping the thin box's left edge
  // (x = 50, y in [5,15]): query box is that degenerate segment's bbox.
  edges.clear();
  index.edgesInBox(50, 8, 50, 12, edges);
  EXPECT_TRUE(std::any_of(edges.begin(), edges.end(), [](const auto& record) {
    return record.obstacle == 2 && record.edge == 3;
  }));

  // Frame corner: a query hugging the hull corner sees the frame's two
  // incident edges and the corner vertex.
  std::vector<SimplifyCandidateIndex::VertexRecord> vertices;
  index.verticesInBox(99, 99, 101, 101, vertices);
  EXPECT_TRUE(std::any_of(vertices.begin(), vertices.end(), [](const auto& record) {
    return record.obstacle == 0 && record.vertex == 2;
  }));
}

TEST(SimplifyCandidateIndex, ApplyRemovalReplacesEdgesAndKeepsIdsStable) {
  // Concave hexagon; remove original vertex 4 (a right-turn corner in ring
  // order); ids never shift.
  const std::vector<Obs> obstacles = {
    makeRing({{0, 0}, {40, 0}, {40, 40}, {20, 40}, {20, 20}, {0, 20}})};
  SimplifyCandidateIndex index(4);
  ASSERT_TRUE(index.build(obstacles));

  // Before: vertex 4 = (20,20) is registered; edges 3 (20,40)->(20,20) and
  // 4 (20,20)->(0,20) exist.
  std::vector<SimplifyCandidateIndex::VertexRecord> vertices;
  index.verticesInBox(20, 20, 20, 20, vertices);
  ASSERT_TRUE(std::any_of(vertices.begin(), vertices.end(), [](const auto& record) {
    return record.obstacle == 0 && record.vertex == 4;
  }));

  index.applyRemoval(0, /*previous=*/3, {20, 40}, /*current=*/4, {20, 20}, /*next=*/5, {0, 20});

  // The removed vertex registration is gone.
  vertices.clear();
  index.verticesInBox(20, 20, 20, 20, vertices);
  EXPECT_FALSE(std::any_of(vertices.begin(), vertices.end(), [](const auto& record) {
    return record.obstacle == 0 && record.vertex == 4;
  }));
  // Edge 4 is gone everywhere it was registered.
  std::vector<SimplifyCandidateIndex::EdgeRecord> edges;
  index.edgesInBox(0, 20, 20, 20, edges);
  EXPECT_FALSE(std::any_of(edges.begin(), edges.end(), [](const auto& record) {
    return record.obstacle == 0 && record.edge == 4;
  }));
  // The chord (20,40)->(0,20) is registered under id 3: its bbox covers
  // (10,30), where neither replaced edge lived.
  edges.clear();
  index.edgesInBox(10, 30, 10, 30, edges);
  EXPECT_TRUE(std::any_of(edges.begin(), edges.end(), [](const auto& record) {
    return record.obstacle == 0 && record.edge == 3;
  }));
  // Dedup cannot prove stale-registration absence (a leftover copy in
  // another bucket would be deduped into the same single report), so pin it
  // geometrically below with a removal whose old-edge bboxes own cells the
  // chord bbox does not touch.
}

TEST(SimplifyCandidateIndex, ApplyRemovalLeavesNoStaleRegistrations) {
  // Pentagon; removing vertex 1 = (6,16) replaces edges 0 and 1 (bboxes
  // reaching y=16) with the horizontal chord (0,10)->(12,10) (bbox y=10
  // only).  Cells above y~12 belonged exclusively to the old edges.
  const std::vector<Obs> obstacles = {makeRing({{0, 10}, {6, 16}, {12, 10}, {12, 0}, {0, 0}})};
  SimplifyCandidateIndex index(4);
  ASSERT_TRUE(index.build(obstacles));

  index.applyRemoval(0, /*previous=*/0, {0, 10}, /*current=*/1, {6, 16}, /*next=*/2, {12, 10});

  // The old edges' exclusive area must hold no registration for ids 0 or 1.
  std::vector<SimplifyCandidateIndex::EdgeRecord> edges;
  index.edgesInBox(0, 14, 12, 16, edges);
  EXPECT_FALSE(std::any_of(edges.begin(), edges.end(), [](const auto& record) {
    return record.obstacle == 0 && (record.edge == 0 || record.edge == 1);
  }));
  // The chord is found on the y=10 strip under the surviving id 0.
  edges.clear();
  index.edgesInBox(5, 10, 7, 10, edges);
  EXPECT_TRUE(std::any_of(edges.begin(), edges.end(), [](const auto& record) {
    return record.obstacle == 0 && record.edge == 0;
  }));
  // Edge 1 (source removed) is gone everywhere.
  edges.clear();
  index.edgesInBox(0, 0, 12, 16, edges);
  EXPECT_FALSE(std::any_of(edges.begin(), edges.end(), [](const auto& record) {
    return record.obstacle == 0 && record.edge == 1;
  }));
}

TEST(SimplifyCandidateIndex, RandomizedRemovalsMatchBruteForceMirror) {
  // Mirrors the index against a live model (alive flags + next links) while
  // applying random removals; after each removal the box differential must
  // still hold with the CURRENT geometry -- this pins chord registration
  // across buckets for arbitrary (non-rectilinear) chords.
  for (std::uint64_t seed = 1; seed <= 30; ++seed) {
    Rng rng{20260926ULL * 7 + seed};
    const auto obstacles = randomRings(rng, rng.range(2, 6));
    SimplifyCandidateIndex index(rng.range(2, 12));
    ASSERT_TRUE(index.build(obstacles));

    // Live mirror in original-id space.
    std::vector<std::vector<char>> alive(obstacles.size());
    std::vector<std::vector<int>> next(obstacles.size());
    std::vector<std::vector<int>> prev(obstacles.size());
    std::vector<int> alive_count(obstacles.size());
    for (size_t o = 0; o < obstacles.size(); ++o) {
      const int n = static_cast<int>(obstacles[o].ordered_vertices_.size());
      alive[o].assign(static_cast<size_t>(n), 1);
      next[o].resize(static_cast<size_t>(n));
      prev[o].resize(static_cast<size_t>(n));
      alive_count[o] = n;
      for (int v = 0; v < n; ++v) {
        next[o][static_cast<size_t>(v)] = (v + 1) % n;
        prev[o][static_cast<size_t>(v)] = (v + n - 1) % n;
      }
    }
    const auto live_brute = [&](int min_x,
                                int min_y,
                                int max_x,
                                int max_y,
                                std::set<std::pair<int, int>>& vertices,
                                std::set<std::pair<int, int>>& edges) {
      for (size_t o = 0; o < obstacles.size(); ++o) {
        const auto& ring = obstacles[o].ordered_vertices_;
        for (size_t v = 0; v < ring.size(); ++v) {
          if (!alive[o][v])
            continue;
          const auto& from = ring[v];
          if (from.first >= min_x && from.first <= max_x && from.second >= min_y &&
              from.second <= max_y)
            vertices.insert({static_cast<int>(o), static_cast<int>(v)});
          const auto& to = ring[static_cast<size_t>(next[o][v])];
          if (std::min(from.first, to.first) <= max_x && min_x <= std::max(from.first, to.first) &&
              std::min(from.second, to.second) <= max_y &&
              min_y <= std::max(from.second, to.second))
            edges.insert({static_cast<int>(o), static_cast<int>(v)});
        }
      }
    };
    const auto check_boxes = [&](int rounds) {
      for (int query = 0; query < rounds; ++query) {
        const int ax = rng.range(-5, 255), ay = rng.range(-5, 255);
        const int bx = rng.range(-5, 255), by = rng.range(-5, 255);
        const int min_x = std::min(ax, bx), max_x = std::max(ax, bx);
        const int min_y = std::min(ay, by), max_y = std::max(ay, by);
        std::set<std::pair<int, int>> expected_vertices, expected_edges;
        live_brute(min_x, min_y, max_x, max_y, expected_vertices, expected_edges);
        std::vector<SimplifyCandidateIndex::VertexRecord> vertices;
        std::vector<SimplifyCandidateIndex::EdgeRecord> edges;
        index.verticesInBox(min_x, min_y, max_x, max_y, vertices);
        index.edgesInBox(min_x, min_y, max_x, max_y, edges);
        std::set<std::pair<int, int>> got_vertices, got_edges;
        for (const auto& record : vertices) got_vertices.insert({record.obstacle, record.vertex});
        for (const auto& record : edges) got_edges.insert({record.obstacle, record.edge});
        for (const auto& expected : expected_vertices)
          ASSERT_TRUE(got_vertices.count(expected)) << "seed " << seed;
        for (const auto& expected : expected_edges)
          ASSERT_TRUE(got_edges.count(expected))
            << "seed " << seed << " edge (" << expected.first << "," << expected.second << ")";
      }
    };

    check_boxes(10);
    for (int removal = 0; removal < 6; ++removal) {
      // Pick a ring that can still lose a vertex, then a random live vertex.
      std::vector<size_t> eligible;
      for (size_t o = 0; o < obstacles.size(); ++o)
        if (alive_count[o] > 3)
          eligible.push_back(o);
      if (eligible.empty())
        break;
      const size_t o =
        eligible[static_cast<size_t>(rng.range(0, static_cast<int>(eligible.size()) - 1))];
      const auto& ring = obstacles[o].ordered_vertices_;
      int b = rng.range(0, static_cast<int>(ring.size()) - 1);
      while (!alive[o][static_cast<size_t>(b)]) b = (b + 1) % static_cast<int>(ring.size());
      const int a = prev[o][static_cast<size_t>(b)];
      const int c = next[o][static_cast<size_t>(b)];
      index.applyRemoval(static_cast<int>(o),
                         a,
                         ring[static_cast<size_t>(a)],
                         b,
                         ring[static_cast<size_t>(b)],
                         c,
                         ring[static_cast<size_t>(c)]);
      alive[o][static_cast<size_t>(b)] = 0;
      next[o][static_cast<size_t>(a)] = c;
      prev[o][static_cast<size_t>(c)] = a;
      --alive_count[o];
      check_boxes(10);
    }
  }
}

TEST(SimplifyCandidateIndex, GateRefusesDuplicateCoordinatesAndZeroLengthEdges) {
  // Duplicate coordinate across rings (Stage 7 unfreeze shape).
  {
    const std::vector<Obs> obstacles = {makeRing({{0, 0}, {10, 0}, {10, 10}, {0, 10}}),
                                        makeRing({{10, 10}, {20, 10}, {20, 20}})};
    SimplifyCandidateIndex index;
    EXPECT_FALSE(index.build(obstacles));
  }
  // Zero-length edge.
  {
    const std::vector<Obs> obstacles = {makeRing({{0, 0}, {10, 0}, {10, 0}, {0, 10}})};
    SimplifyCandidateIndex index;
    EXPECT_FALSE(index.build(obstacles));
  }
  // Tombstoned empty rings are fine.
  {
    std::vector<Obs> obstacles = {makeRing({}), makeRing({{0, 0}, {10, 0}, {10, 10}, {0, 10}})};
    SimplifyCandidateIndex index;
    EXPECT_TRUE(index.build(obstacles));
  }
}

}  // namespace
}  // namespace raystar
