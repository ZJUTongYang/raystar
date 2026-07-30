#include <gtest/gtest.h>
#include <raystar/raystar_core.h>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

using namespace raystar;

static GridMap makeTestmapFromPgm() {
  GridMap map;
  map.width = 50;
  map.height = 50;
  map.resolution = 0.1f;
  map.origin_x = 2.0;
  map.origin_y = 3.0;
  map.data.resize(2500, 0);

  FILE* f = fopen(RAYSTAR_TESTMAP_PATH, "rb");
  if (!f) {
    map.width = 0;
    return map;
  }
  fseek(f, 0, SEEK_END);
  long sz = ftell(f);
  fseek(f, 0, SEEK_SET);
  std::vector<unsigned char> raw(sz);
  const size_t bytes_read = fread(raw.data(), 1, static_cast<size_t>(sz), f);
  fclose(f);
  if (bytes_read != static_cast<size_t>(sz)) {
    map.width = 0;
    map.height = 0;
    map.data.clear();
    return map;
  }

  size_t pos = 0;
  int nl = 0;
  while (pos < raw.size() && nl < 3) {
    if (raw[pos] == '\n')
      nl++;
    pos++;
  }
  for (size_t i = 0; i < 2500 && pos < raw.size(); ++i, ++pos) {
    map.data[i] = (raw[pos] > 200) ? 0 : 1;
  }
  return map;
}

static const BoundaryEndpoint* findEndpoint(const VisibilityRegion& region,
                                            double x,
                                            double y,
                                            double tolerance = 1e-6) {
  const auto it = std::find_if(region.begin(), region.end(), [&](const auto& endpoint) {
    return std::abs(endpoint.position.first - x) < tolerance &&
           std::abs(endpoint.position.second - y) < tolerance;
  });
  return it == region.end() ? nullptr : &*it;
}

static Polymap makeReadyPolymap(
  const GridMap& map, int start_x, int start_y, int goal_x, int goal_y) {
  auto created = Polymap::create(map, start_x, start_y, goal_x, goal_y);
  if (!created) {
    throw std::runtime_error("Polymap construction failed with status " +
                             std::to_string(static_cast<int>(created.status)) + ": " +
                             created.error);
  }
  return std::move(*created.value);
}

TEST(TopoVRegression, LocateVertexCorrect) {
  GridMap map = makeTestmapFromPgm();
  ASSERT_EQ(map.width, 50u);

  Polymap pm = makeReadyPolymap(map, 2, 2, 46, 46);
  ASSERT_TRUE(pm.hasSolution());

  // The first obstacle discovered by getPolyObstacles (BFS from start) is the
  // L-shaped block at x=9..14, y=7..15. After simplifyPolyObstacles it has 5
  // vertices: obs[0] = (9,16)(9,7)(10,7)(13,14)(13,16).
  auto v1 = pm.locateVertex(10, 7);
  EXPECT_EQ(v1.first, 0);
  EXPECT_EQ(v1.second, 2) << "locateVertex(10,7) should return (0,2) since obs[0][2]=(10,7)";

  auto v2 = pm.locateVertex(9, 7);
  EXPECT_EQ(v2.first, 0);
  EXPECT_EQ(v2.second, 1) << "locateVertex(9,7) should return (0,1) since obs[0][1]=(9,7)";
}

TEST(TopoVRegression, VisibilityEndpointMetadataMatchesGeometry) {
  GridMap map = makeTestmapFromPgm();
  map.data[2 * map.width + 2] = 0;
  for (unsigned int x = 0; x < map.width; ++x) {
    map.data[x] = 1;
    map.data[(map.height - 1) * map.width + x] = 1;
  }
  for (unsigned int y = 0; y < map.height; ++y) {
    map.data[y * map.width] = 1;
    map.data[y * map.width + map.width - 1] = 1;
  }
  Polymap pm = makeReadyPolymap(map, 2, 2, 46, 46);
  ASSERT_TRUE(pm.hasSolution());

  const auto expect_edge = [&](const BoundaryEndpoint* endpoint,
                               std::pair<int, int> expected_from,
                               std::pair<int, int> expected_to) {
    ASSERT_NE(endpoint, nullptr);
    const auto* edge = std::get_if<DirectedObstacleEdge>(&endpoint->support);
    ASSERT_NE(edge, nullptr) << "endpoint must retain supporting-edge metadata";
    ASSERT_TRUE(pm.isValidTopology({edge->from.obstacle, edge->from.vertex}));
    ASSERT_TRUE(pm.isValidTopology({edge->to.obstacle, edge->to.vertex}));
    const auto& from = pm.obstacles()[static_cast<size_t>(edge->from.obstacle)]
                         .ordered_vertices_[static_cast<size_t>(edge->from.vertex)];
    const auto& to = pm.obstacles()[static_cast<size_t>(edge->to.obstacle)]
                       .ordered_vertices_[static_cast<size_t>(edge->to.vertex)];
    EXPECT_EQ(from, expected_from);
    EXPECT_EQ(to, expected_to);
  };

  VisibilityRegion from_root;
  std::string error;
  ASSERT_TRUE(pm.getVisibilityRegion(2, 2, from_root, &error)) << error;
  for (const auto& endpoint : from_root) {
    EXPECT_TRUE(pm.validateBoundaryEndpoint(endpoint, &error)) << error;
    EXPECT_FALSE(std::holds_alternative<std::monostate>(endpoint.support));
    EXPECT_TRUE(endpoint.exact_position.has_value());
  }

  expect_edge(findEndpoint(from_root, 14.199273, 41.924894, 1e-5), {1, 40}, {49, 47});
  const auto* exact_rational_intersection = findEndpoint(from_root, 15.0, 10.125);
  expect_edge(exact_rational_intersection, {15, 12}, {15, 8});
  ASSERT_NE(exact_rational_intersection, nullptr);
  EXPECT_EQ(exactPoint(*exact_rational_intersection).x(), exact_geometry::FT(15));
  EXPECT_EQ(exactPoint(*exact_rational_intersection).y() * exact_geometry::FT(8),
            exact_geometry::FT(81));

  VisibilityRegion from_obstacle_vertex;
  ASSERT_TRUE(pm.getVisibilityRegion(10, 7, from_obstacle_vertex, &error)) << error;
  for (const auto& endpoint : from_obstacle_vertex) {
    EXPECT_TRUE(pm.validateBoundaryEndpoint(endpoint, &error)) << error;
    EXPECT_FALSE(std::holds_alternative<std::monostate>(endpoint.support));
    EXPECT_TRUE(endpoint.exact_position.has_value());
  }

  const auto* integer_edge_point = findEndpoint(from_obstacle_vertex, 17.0, 14.0);
  expect_edge(integer_edge_point, {15, 14}, {18, 14});
  EXPECT_EQ(pm.locateVertex(17.0, 14.0), std::make_pair(-1, -1));
  ASSERT_NE(integer_edge_point, nullptr);
  const auto* horizontal_edge = std::get_if<DirectedObstacleEdge>(&integer_edge_point->support);
  ASSERT_NE(horizontal_edge, nullptr);
  const BoundaryEndpoint earlier_on_edge{{16.0, 14.0}, *horizontal_edge};
  EXPECT_TRUE(pm.boundarySupportsConsecutive(earlier_on_edge, *integer_edge_point));
  EXPECT_FALSE(pm.boundarySupportsConsecutive(*integer_edge_point, earlier_on_edge));
  expect_edge(findEndpoint(from_obstacle_vertex, 30.0, 9.857142857, 1e-6), {30, 16}, {30, 7});

  std::vector<std::pair<double, double>> legacy_positions;
  std::vector<std::pair<int, int>> legacy_topology;
  ASSERT_TRUE(pm.getVisibilityRegion(10, 7, legacy_positions, legacy_topology, &error)) << error;
  ASSERT_EQ(legacy_positions.size(), from_obstacle_vertex.size());
  ASSERT_EQ(legacy_topology.size(), from_obstacle_vertex.size());
  const auto integer_edge_position =
    std::find(legacy_positions.begin(), legacy_positions.end(), std::make_pair(17.0, 14.0));
  ASSERT_NE(integer_edge_position, legacy_positions.end());
  const size_t integer_edge_index =
    static_cast<size_t>(std::distance(legacy_positions.begin(), integer_edge_position));
  EXPECT_EQ(legacy_topology[integer_edge_index], std::make_pair(-1, -1));
}

TEST(TopoVRegression, Node0TopoVMatchesLocateVertex) {
  GridMap map = makeTestmapFromPgm();
  RaystarCore core;
  auto result = core.plan(map, 2, 2, 46, 46, 3, false);
  ASSERT_TRUE(result.success);

  const auto& nodes = core.getNodes();
  ASSERT_GE(nodes.size(), 2u);

  const auto& pm = *result.polymap;
  ASSERT_EQ(nodes[0].visibility().size(), nodes[0].visibilityTopology().size());
  int mismatch_count = 0;
  for (size_t i = 0; i < nodes[0].visibility().size(); ++i) {
    int ix = static_cast<int>(nodes[0].visibility()[i].first);
    int iy = static_cast<int>(nodes[0].visibility()[i].second);
    bool is_integer =
      (nodes[0].visibility()[i].first == ix) && (nodes[0].visibility()[i].second == iy);
    if (!is_integer || nodes[0].visibilityTopology()[i].first < 0)
      continue;

    auto lv = pm.locateVertex(ix, iy);
    if (lv.first != nodes[0].visibilityTopology()[i].first ||
        lv.second != nodes[0].visibilityTopology()[i].second) {
      ADD_FAILURE() << "V_[" << i << "]=(" << ix << "," << iy << ") topo=("
                    << nodes[0].visibilityTopology()[i].first << ","
                    << nodes[0].visibilityTopology()[i].second << ") but locateVertex=(" << lv.first
                    << "," << lv.second << ")";
      mismatch_count++;
    }
  }
  EXPECT_EQ(mismatch_count, 0)
    << "All integer obstacle vertices in Node 0's V_ should have matching topo_V_";
}

TEST(TopoVRegression, Node1EndAngleReasonable) {
  GridMap map = makeTestmapFromPgm();
  RaystarCore core;
  auto result = core.plan(map, 2, 2, 46, 46, 3, false);
  ASSERT_TRUE(result.success);

  const auto& nodes = core.getNodes();
  ASSERT_GE(nodes.size(), 2u);

  double angle_span = nodes[1].endAngle() - nodes[1].startAngle();
  if (angle_span < 0)
    angle_span += kTwoPi;

  EXPECT_LT(angle_span, kPi) << "Node 1 angle span should be < pi (180deg) for a gap, got "
                             << angle_span << " rad (" << angle_span * 180 / kPi << " deg)";
}

TEST(TopoVRegression, Node1FullVisibilityGeometry) {
  GridMap map = makeTestmapFromPgm();
  RaystarCore core;
  auto result = core.plan(map, 2, 2, 46, 46, 3, false);
  ASSERT_TRUE(result.success);

  const auto& nodes = core.getNodes();
  ASSERT_GE(nodes.size(), 2u);
  ASSERT_EQ(nodes[1].seed(), std::make_pair(10, 7));

  const auto& full_v = nodes[1].fullVisibility();
  ASSERT_GE(full_v.size(), 16u);

  const std::vector<std::pair<double, double>> expected_prefix = {
    {9.0, 7.0}, {1.0, 7.0}, {1.0, 2.0}, {40.0, 1.0}, {41.17391304347826, 7.0}, {30.0, 7.0}};
  for (size_t i = 0; i < expected_prefix.size(); ++i) {
    EXPECT_NEAR(full_v[i].first, expected_prefix[i].first, 1e-6) << "vertex " << i;
    EXPECT_NEAR(full_v[i].second, expected_prefix[i].second, 1e-6) << "vertex " << i;
  }

  auto count_near = [&](double x, double y) {
    return std::count_if(full_v.begin(), full_v.end(), [&](const auto& v) {
      return std::abs(v.first - x) < 1e-6 && std::abs(v.second - y) < 1e-6;
    });
  };
  EXPECT_EQ(count_near(9.0, 16.0), 0)
    << "A vertex inside the source obstacle's blocked wedge must not be emitted";
  EXPECT_EQ(count_near(13.0, 14.0), 1)
    << "The open visibility fan's closing vertex must be emitted exactly once";

  const exact_geometry::Point source(10, 7);
  const exact_geometry::Point start(9, 7);
  const exact_geometry::Point end(13, 14);
  for (const auto& endpoint : nodes[1].fullVisibilityRegion()) {
    EXPECT_TRUE(
      exact_geometry::isClosedCounterClockwiseSweepMember(source, start, end, exactPoint(endpoint)))
      << "Visibility vertex (" << endpoint.position.first << "," << endpoint.position.second
      << ") lies inside the source obstacle's blocked wedge";
  }
}

TEST(TopoVRegression, ScopedVisibilityPreservesBoundarySupport) {
  GridMap map = makeTestmapFromPgm();
  RaystarCore core;
  auto result = core.plan(map, 2, 2, 46, 46, 3, false);
  ASSERT_TRUE(result.success) << result.message;

  const Node* node = nullptr;
  for (const auto& candidate : core.getNodes()) {
    if (candidate.seed() == std::make_pair(10, 7)) {
      node = &candidate;
      break;
    }
  }
  ASSERT_NE(node, nullptr);
  ASSERT_EQ(node->parentIndex(), 0);

  const auto* scoped_fractional = findEndpoint(node->visibilityRegion(), 15.0, 10.125);
  const auto* parent_fractional =
    findEndpoint(core.getNodes().front().visibilityRegion(), 15.0, 10.125);
  ASSERT_NE(scoped_fractional, nullptr);
  ASSERT_NE(parent_fractional, nullptr);
  ASSERT_NE(std::get_if<DirectedObstacleEdge>(&scoped_fractional->support), nullptr);
  EXPECT_TRUE(scoped_fractional->support == parent_fractional->support)
    << "scoping must copy the complete fractional endpoint";

  const auto* closing_vertex = findEndpoint(node->visibilityRegion(), 13.0, 14.0);
  ASSERT_NE(closing_vertex, nullptr);
  const auto* exact = std::get_if<ObstacleVertexId>(&closing_vertex->support);
  ASSERT_NE(exact, nullptr);
  EXPECT_EQ(result.polymap->locateVertex(13, 14), std::make_pair(exact->obstacle, exact->vertex));
}

TEST(TopoVRegression, ExactSameRayOrderingAvoidsSpuriousWideScope) {
  GridMap map = makeTestmapFromPgm();
  RaystarCore core;
  auto result = core.plan(map, 2, 2, 46, 46, 3, false);
  ASSERT_TRUE(result.success);

  const auto& nodes = core.getNodes();
  const Node* node = nullptr;
  size_t matching_node_count = 0;
  for (const auto& candidate : nodes) {
    if (candidate.seed() != std::make_pair(15, 16) || candidate.parentIndex() < 0 ||
        candidate.parentIndex() >= static_cast<int>(nodes.size()) ||
        nodes[static_cast<size_t>(candidate.parentIndex())].seed() != std::make_pair(10, 7)) {
      continue;
    }
    node = &candidate;
    ++matching_node_count;
  }
  ASSERT_EQ(matching_node_count, 1u);
  ASSERT_NE(node, nullptr);

  const auto* near = findEndpoint(node->visibilityRegion(), 24.0, 23.0);
  const auto* far = findEndpoint(node->visibilityRegion(), 26.57142857142857, 25.0);
  ASSERT_NE(near, nullptr);
  ASSERT_NE(far, nullptr);
  const exact_geometry::Point source(15, 16);
  EXPECT_TRUE(exact_geometry::isSameDirectedRay(source, exactPoint(*near), exactPoint(*far)));

  const auto near_index = static_cast<size_t>(near - node->visibilityRegion().data());
  const auto far_index = static_cast<size_t>(far - node->visibilityRegion().data());
  EXPECT_LT(near_index, far_index)
    << "Exact equal-ray sorting must retain the boundary-list near/far order";

  const auto child =
    std::find_if(node->children().begin(), node->children().end(), [&](const Child& candidate) {
      return exactPoint(candidate.endpoint()) == exactPoint(*near) &&
             exactPoint(candidate.oppositeEndpoint()) == exactPoint(*far);
    });
  ASSERT_NE(child, node->children().end());
  EXPECT_FALSE(child->isLeftGap())
    << "A floating atan2 perturbation must not turn this right gap into a wide left scope";
}

TEST(TopoVRegression, SameRayBoundaryOrderPreserved) {
  GridMap map = makeTestmapFromPgm();
  RaystarCore core;
  auto result = core.plan(map, 2, 2, 46, 46, 3, false);
  ASSERT_TRUE(result.success);

  const Node* node = nullptr;
  for (const auto& candidate : core.getNodes()) {
    if (candidate.seed() == std::make_pair(15, 12)) {
      node = &candidate;
      break;
    }
  }
  ASSERT_NE(node, nullptr);
  ASSERT_GE(node->fullVisibility().size(), 2u);

  EXPECT_NEAR(node->fullVisibility()[0].first, 17.0, 1e-6);
  EXPECT_NEAR(node->fullVisibility()[0].second, 12.0, 1e-6);
  EXPECT_NEAR(node->fullVisibility()[1].first, 30.0, 1e-6);
  EXPECT_NEAR(node->fullVisibility()[1].second, 12.0, 1e-6);
}

TEST(TopoVRegression, ChildCtopoMatchesLocateVertex) {
  GridMap map = makeTestmapFromPgm();
  RaystarCore core;
  auto result = core.plan(map, 2, 2, 46, 46, 3, false);
  ASSERT_TRUE(result.success);

  const auto& nodes = core.getNodes();
  ASSERT_GE(nodes.size(), 1u);
  ASSERT_GT(nodes[0].children().size(), 0u);

  const auto& pm = *result.polymap;
  bool found_fractional_far_endpoint = false;
  for (size_t ni = 0; ni < nodes.size(); ++ni) {
    for (size_t ci = 0; ci < nodes[ni].children().size(); ++ci) {
      const auto& c = nodes[ni].children()[ci];
      std::string endpoint_error;
      EXPECT_TRUE(pm.validateBoundaryEndpoint(c.endpoint(), &endpoint_error)) << endpoint_error;
      EXPECT_TRUE(pm.validateBoundaryEndpoint(c.oppositeEndpoint(), &endpoint_error))
        << endpoint_error;
      const auto* exact = std::get_if<ObstacleVertexId>(&c.endpoint().support);
      ASSERT_NE(exact, nullptr) << "node " << ni << " child " << ci;
      auto lv = pm.locateVertex(c.coordinate().first, c.coordinate().second);
      EXPECT_EQ(lv.first, c.obstacleIndex())
        << "node " << ni << " child " << ci << " c_=(" << c.coordinate().first << ","
        << c.coordinate().second << ") c_obs_index mismatch";
      EXPECT_EQ(lv.second, c.vertexIndex())
        << "node " << ni << " child " << ci << " c_=(" << c.coordinate().first << ","
        << c.coordinate().second << ") c_ver_index mismatch";
      EXPECT_EQ(lv, std::make_pair(exact->obstacle, exact->vertex));
      EXPECT_EQ(c.oppositeCoordinate(), c.oppositeEndpoint().position);
      EXPECT_FALSE(std::holds_alternative<std::monostate>(c.oppositeEndpoint().support));

      const bool is_fractional = std::abs(c.oppositeEndpoint().position.first -
                                          std::round(c.oppositeEndpoint().position.first)) > 1e-9 ||
                                 std::abs(c.oppositeEndpoint().position.second -
                                          std::round(c.oppositeEndpoint().position.second)) > 1e-9;
      if (is_fractional &&
          std::holds_alternative<DirectedObstacleEdge>(c.oppositeEndpoint().support)) {
        found_fractional_far_endpoint = true;
        EXPECT_EQ(c.oppositeObstacleIndex(), -1);
        EXPECT_EQ(c.oppositeVertexIndex(), -1);
      }
    }
  }
  EXPECT_TRUE(found_fractional_far_endpoint)
    << "At least one generated child must retain a fractional far-endpoint edge support";
}
