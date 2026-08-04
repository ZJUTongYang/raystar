#include <gtest/gtest.h>
#include <raystar/raystar_core.h>
#include <exact_point_location.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <initializer_list>
#include <limits>
#include <mutex>
#include <thread>
#include <type_traits>

using namespace raystar;

static BoundaryEndpoint makeExactBoundaryEndpoint(const exact_geometry::Point& point) {
  BoundaryEndpoint endpoint;
  endpoint.position = {CGAL::to_double(point.x()), CGAL::to_double(point.y())};
  endpoint.exact_position = point;
  endpoint.support = std::monostate{};
  return endpoint;
}

static VisibilityRegion makeIntegerBoundary(std::initializer_list<std::pair<int, int>> vertices) {
  VisibilityRegion boundary;
  boundary.reserve(vertices.size());
  for (const auto& vertex : vertices) {
    boundary.emplace_back(
      makeExactBoundaryEndpoint(exact_geometry::Point(vertex.first, vertex.second)));
  }
  return boundary;
}

TEST(ExactPointLocation, ClassifiesClosedRectangleAndBoundary) {
  const VisibilityRegion rectangle = makeIntegerBoundary({{0, 0}, {4, 0}, {4, 4}, {0, 4}});
  const auto classify = [&](int x, int y) {
    return detail::classifyPointInVisibilityBoundary(rectangle, exact_geometry::Point(x, y));
  };

  EXPECT_EQ(classify(2, 2), std::make_pair(true, false));
  EXPECT_EQ(classify(-1, 2), std::make_pair(false, false));
  EXPECT_EQ(classify(5, 2), std::make_pair(false, false));
  EXPECT_EQ(classify(0, 2), std::make_pair(false, true));
  EXPECT_EQ(classify(4, 2), std::make_pair(false, true));
  EXPECT_EQ(classify(2, 0), std::make_pair(false, true));
  EXPECT_EQ(classify(2, 4), std::make_pair(false, true));
  EXPECT_EQ(classify(0, 0), std::make_pair(false, true));
  EXPECT_EQ(classify(4, 4), std::make_pair(false, true));

  VisibilityRegion reversed = rectangle;
  std::reverse(reversed.begin(), reversed.end());
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(reversed, exact_geometry::Point(2, 2)),
            std::make_pair(true, false));
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(reversed, exact_geometry::Point(5, 2)),
            std::make_pair(false, false));
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(VisibilityRegion{rectangle[0], rectangle[1]},
                                                      exact_geometry::Point(2, 0)),
            std::make_pair(false, false));
}

TEST(ExactPointLocation, CooperativeOverloadDistinguishesFailureAndStop) {
  const VisibilityRegion rectangle = makeIntegerBoundary({{0, 0}, {4, 0}, {4, 4}, {0, 4}});
  const StopToken never_stop;
  std::pair<bool, bool> classification{false, true};

  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(
              rectangle, exact_geometry::Point(2, 2), classification, never_stop),
            OperationStatus::success);
  EXPECT_EQ(classification, std::make_pair(true, false));

  const VisibilityRegion invalid_boundary = makeIntegerBoundary({{0, 0}, {4, 0}});
  classification = {true, true};
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(
              invalid_boundary, exact_geometry::Point(2, 1), classification, never_stop),
            OperationStatus::failure);
  EXPECT_EQ(classification, std::make_pair(true, true));

  VisibilityRegion long_boundary;
  long_boundary.reserve(256);
  for (int x = 0; x < 256; ++x) {
    long_boundary.emplace_back(makeExactBoundaryEndpoint(exact_geometry::Point(x, 0)));
  }

  size_t poll_count = 0;
  const StopToken stop_token([&poll_count]() { return ++poll_count == 8; });
  classification = {true, true};
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(
              long_boundary, exact_geometry::Point(0, 1), classification, stop_token),
            OperationStatus::stopped);
  EXPECT_EQ(classification, std::make_pair(true, true));
  EXPECT_EQ(poll_count, 8u);
}

TEST(ExactPointLocation, ExplicitClosedModePreservesRectangleSemantics) {
  const VisibilityRegion rectangle = makeIntegerBoundary({{0, 0}, {4, 0}, {4, 4}, {0, 4}});
  const exact_geometry::Point unused_source(10, 10);
  const auto classify = [&](int x, int y) {
    const exact_geometry::Point query(x, y);
    const auto compatibility_result = detail::classifyPointInVisibilityBoundary(rectangle, query);
    const auto explicit_result = detail::classifyPointInVisibilityBoundary(
      rectangle, query, unused_source, detail::VisibilityBoundaryMode::closed_cycle);
    EXPECT_EQ(explicit_result, compatibility_result);
    return explicit_result;
  };

  EXPECT_EQ(classify(2, 2), std::make_pair(true, false));
  EXPECT_EQ(classify(5, 2), std::make_pair(false, false));
  EXPECT_EQ(classify(0, 2), std::make_pair(false, true));
  EXPECT_EQ(classify(0, 0), std::make_pair(false, true));
}

TEST(ExactPointLocation, OpenSectorClosesThroughSourceSpokes) {
  const exact_geometry::Point source(0, 0);
  const VisibilityRegion open_boundary = makeIntegerBoundary({{0, 4}, {4, 4}, {4, 0}});

  // Closing through source gives the square [0, 4] x [0, 4].  The old
  // artificial last-to-first chord instead forms the upper-right triangle
  // and incorrectly excludes this point near source.
  const exact_geometry::Point near_source(1, 1);
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(
              open_boundary, near_source, source, detail::VisibilityBoundaryMode::open_sector),
            std::make_pair(true, false));
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(open_boundary, near_source),
            std::make_pair(false, false));

  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(
              open_boundary, source, source, detail::VisibilityBoundaryMode::open_sector),
            std::make_pair(false, true));
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(open_boundary,
                                                      exact_geometry::Point(0, 2),
                                                      source,
                                                      detail::VisibilityBoundaryMode::open_sector),
            std::make_pair(false, true));
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(open_boundary,
                                                      exact_geometry::Point(2, 0),
                                                      source,
                                                      detail::VisibilityBoundaryMode::open_sector),
            std::make_pair(false, true));
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(open_boundary,
                                                      exact_geometry::Point(5, 2),
                                                      source,
                                                      detail::VisibilityBoundaryMode::open_sector),
            std::make_pair(false, false));
}

TEST(ExactPointLocation, UsesExactEndpointsAtLargeCoordinates) {
  using exact_geometry::FT;
  using exact_geometry::Point;
  const FT coordinate(2000000000);
  const FT delta = FT(1) / FT(1000000000);

  const auto make_rectangle = [&](const FT& right_x) {
    return VisibilityRegion{makeExactBoundaryEndpoint(Point(coordinate - FT(2), FT(0))),
                            makeExactBoundaryEndpoint(Point(right_x, FT(0))),
                            makeExactBoundaryEndpoint(Point(right_x, FT(2))),
                            makeExactBoundaryEndpoint(Point(coordinate - FT(2), FT(2)))};
  };
  const VisibilityRegion containing = make_rectangle(coordinate + delta);
  const VisibilityRegion excluding = make_rectangle(coordinate - delta);

  ASSERT_EQ(containing[1].position.first, excluding[1].position.first)
    << "The projected doubles intentionally lose the exact side of the query";
  const Point query(2000000000, 1);
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(containing, query),
            std::make_pair(true, false));
  EXPECT_EQ(detail::classifyPointInVisibilityBoundary(excluding, query),
            std::make_pair(false, false));
}

TEST(ExactPointLocation, HandlesWeaklySimpleSameRayBacktrack) {
  VisibilityRegion boundary =
    makeIntegerBoundary({{0, 0}, {6, 0}, {6, 6}, {3, 6}, {3, 3}, {3, 6}, {0, 6}});
  const auto expect_classification = [&](const VisibilityRegion& candidate) {
    EXPECT_EQ(detail::classifyPointInVisibilityBoundary(candidate, exact_geometry::Point(1, 3)),
              std::make_pair(true, false));
    EXPECT_EQ(detail::classifyPointInVisibilityBoundary(candidate, exact_geometry::Point(4, 4)),
              std::make_pair(true, false));
    EXPECT_EQ(detail::classifyPointInVisibilityBoundary(candidate, exact_geometry::Point(7, 4)),
              std::make_pair(false, false));
    EXPECT_EQ(detail::classifyPointInVisibilityBoundary(candidate, exact_geometry::Point(3, 4)),
              std::make_pair(false, true));
  };

  expect_classification(boundary);
  std::reverse(boundary.begin(), boundary.end());
  expect_classification(boundary);
}

TEST(ExactPointLocation, MatchesLegacyParityAwayFromBoundary) {
  const std::vector<VisibilityRegion> boundaries = {
    makeIntegerBoundary({{0, 0}, {8, 0}, {8, 8}, {0, 8}}),
    makeIntegerBoundary({{0, 0}, {8, 0}, {4, 0}, {8, 4}, {8, 8}, {0, 8}}),
    makeIntegerBoundary({{0, 0}, {8, 0}, {8, 0}, {8, 8}, {4, 8}, {4, 4}, {4, 8}, {0, 8}})};

  const auto legacy_classify = [](
                                 const VisibilityRegion& boundary, double query_x, double query_y) {
    bool right_ray = false;
    bool left_ray = false;
    for (size_t i = 0, j = boundary.size() - 1; i < boundary.size(); j = i++) {
      const auto& current = boundary[i].position;
      const auto& previous = boundary[j].position;
      const bool crosses = (current.second > query_y) != (previous.second > query_y);
      if (!crosses)
        continue;
      const double intersection_x = current.first + (previous.first - current.first) *
                                                      (query_y - current.second) /
                                                      (previous.second - current.second);
      if (query_x < intersection_x)
        right_ray = !right_ray;
      if (query_x > intersection_x)
        left_ray = !left_ray;
    }
    if (right_ray == left_ray)
      return std::make_pair(right_ray, false);
    return std::make_pair(false, true);
  };

  for (size_t boundary_index = 0; boundary_index < boundaries.size(); ++boundary_index) {
    for (int x_quarters = -3; x_quarters <= 35; x_quarters += 2) {
      for (int y_quarters = -3; y_quarters <= 35; y_quarters += 2) {
        SCOPED_TRACE(testing::Message() << "boundary=" << boundary_index << " query_quarters=("
                                        << x_quarters << "," << y_quarters << ")");
        const exact_geometry::Point exact_query(
          exact_geometry::FT(x_quarters) / exact_geometry::FT(4),
          exact_geometry::FT(y_quarters) / exact_geometry::FT(4));
        const auto exact_result =
          detail::classifyPointInVisibilityBoundary(boundaries[boundary_index], exact_query);
        if (exact_result.second)
          continue;
        EXPECT_EQ(exact_result,
                  legacy_classify(boundaries[boundary_index], x_quarters / 4.0, y_quarters / 4.0));
      }
    }
  }
}

static GridMap makeOpenMap() {
  GridMap map;
  map.width = 20;
  map.height = 20;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.resize(400, 0);
  return map;
}

static GridMap makeSimpleMap() {
  GridMap map;
  map.width = 30;
  map.height = 30;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.resize(900, 0);

  for (unsigned int y = 10; y < 20; ++y)
    for (unsigned int x = 10; x < 20; ++x) map.data[y * 30 + x] = 1;

  return map;
}

TEST(RaystarCore, PlanStraightPath) {
  auto map = makeOpenMap();

  RaystarCore core;
  auto result = core.plan(map, 2, 2, 17, 17, 1, false);

  EXPECT_TRUE(result.success);
  EXPECT_EQ(result.outcome, PlanningOutcome::complete);
  EXPECT_EQ(result.completion, PlanningCompletion::requested_k_reached);
  EXPECT_GE(result.path_solutions.size(), 1u);
  EXPECT_EQ(result.expanded_nodes, core.getNodes().size());
  EXPECT_GT(result.path_solutions[0].path_cost_, 0.0);
  for (const auto& node : core.getNodes()) {
    EXPECT_GE(node.V_.size(), 2u);
    EXPECT_EQ(node.V_.size(), node.topo_V_.size());
  }
}

TEST(RaystarCore, EnumeratesEveryPathWithinInclusiveCostBound) {
  const auto map = makeSimpleMap();
  RaystarCore core;
  const auto baseline = core.plan(map, 5, 14, 25, 15, 2, false);
  ASSERT_TRUE(baseline.success) << baseline.message;
  ASSERT_EQ(baseline.path_solutions.size(), 2u);
  ASSERT_LT(baseline.path_solutions[0].path_cost_, baseline.path_solutions[1].path_cost_);

  const double between =
    (baseline.path_solutions[0].path_cost_ + baseline.path_solutions[1].path_cost_) / 2.0;
  const auto first_only =
    core.plan(map, 5, 14, 25, 15, SearchObjective::allWithinCost(between), false);
  ASSERT_TRUE(first_only.success) << first_only.message;
  ASSERT_EQ(first_only.path_solutions.size(), 1u);
  EXPECT_DOUBLE_EQ(first_only.path_solutions.front().path_cost_,
                   baseline.path_solutions.front().path_cost_);
  EXPECT_EQ(first_only.completion, PlanningCompletion::cost_bound_exhausted);
  EXPECT_EQ(first_only.limit_reached, PlanningLimitReached::none);

  const double inclusive_bound = baseline.path_solutions.back().path_cost_;
  const auto both =
    core.plan(map, 5, 14, 25, 15, SearchObjective::allWithinCost(inclusive_bound), false);
  ASSERT_TRUE(both.success) << both.message;
  ASSERT_EQ(both.path_solutions.size(), baseline.path_solutions.size());
  for (size_t i = 0; i < both.path_solutions.size(); ++i) {
    EXPECT_DOUBLE_EQ(both.path_solutions[i].path_cost_, baseline.path_solutions[i].path_cost_);
    EXPECT_EQ(both.path_solutions[i].projectedPath(), baseline.path_solutions[i].projectedPath());
    EXPECT_LE(both.path_solutions[i].path_cost_, inclusive_bound);
  }
  EXPECT_TRUE(both.completion == PlanningCompletion::cost_bound_exhausted ||
              both.completion == PlanningCompletion::frontier_exhausted);
}

TEST(RaystarCore, CostBoundBelowShortestPathCompletesWithNoResult) {
  const auto map = makeOpenMap();
  RaystarCore core;
  const double direct_cost = std::hypot(15.0, 15.0);

  const auto below = core.plan(
    map, 2, 2, 17, 17, SearchObjective::allWithinCost(std::nextafter(direct_cost, 0.0)), false);
  EXPECT_FALSE(below.success);
  EXPECT_EQ(below.outcome, PlanningOutcome::no_path);
  // A cheap certified frontier lower bound may overlap this one-ULP request,
  // requiring expansion of the root before exact solution comparison rejects
  // the path. Both proof routes are complete.
  EXPECT_TRUE(below.completion == PlanningCompletion::cost_bound_exhausted ||
              below.completion == PlanningCompletion::frontier_exhausted);
  EXPECT_TRUE(below.path_solutions.empty());
  EXPECT_LE(core.getNodes().size(), 1u);

  const auto equal =
    core.plan(map, 2, 2, 17, 17, SearchObjective::allWithinCost(direct_cost), false);
  ASSERT_TRUE(equal.success) << equal.message;
  ASSERT_EQ(equal.path_solutions.size(), 1u);
  EXPECT_DOUBLE_EQ(equal.path_solutions.front().path_cost_, direct_cost);
  EXPECT_TRUE(equal.completion == PlanningCompletion::cost_bound_exhausted ||
              equal.completion == PlanningCompletion::frontier_exhausted);
}

TEST(RaystarCore, RadicalSumBoundaryIsCompleteForSingleAndSharedTreeSearch) {
  GridMap map;
  map.width = 60;
  map.height = 60;
  map.resolution = 1.0f;
  map.data.assign(60U * 60U, 0);
  map.data[30U * 60U + 20U] = 1;

  const PlanEndpoint start(19, 52, 19.5, 52.5);
  const PlanEndpoint goal(20, 29, 20.5, 29.0);
  constexpr double inclusive_bound = 0x1.79fa384f9da53p+4;
  RaystarCore core;

  const auto single =
    core.plan(map, start, goal, SearchObjective::allWithinCost(inclusive_bound), false);
  ASSERT_TRUE(single.success) << single.message;
  ASSERT_EQ(single.path_solutions.size(), 1u);
  EXPECT_DOUBLE_EQ(single.path_solutions.front().path_cost_, inclusive_bound);
  ASSERT_EQ(single.path_solutions.front().turning_points_.size(), 1u);
  EXPECT_EQ(single.path_solutions.front().turning_points_.front(), std::make_pair(20, 30));
  EXPECT_TRUE(single.completion == PlanningCompletion::cost_bound_exhausted ||
              single.completion == PlanningCompletion::frontier_exhausted);

  const auto multi = core.planToGoalsWithinCosts(map, start, {{goal, inclusive_bound}}, false);
  ASSERT_EQ(multi.goal_results.size(), 1u);
  const auto& goal_result = multi.goal_results.front();
  ASSERT_TRUE(goal_result.success) << goal_result.message;
  ASSERT_EQ(goal_result.path_solutions.size(), 1u);
  EXPECT_DOUBLE_EQ(goal_result.path_solutions.front().path_cost_, inclusive_bound);
  EXPECT_EQ(goal_result.path_solutions.front().projectedPath(),
            single.path_solutions.front().projectedPath());

  const double excluded_bound = std::nextafter(inclusive_bound, 0.0);
  const auto excluded =
    core.plan(map, start, goal, SearchObjective::allWithinCost(excluded_bound), false);
  EXPECT_FALSE(excluded.success);
  EXPECT_TRUE(excluded.path_solutions.empty());
  EXPECT_EQ(excluded.completion, PlanningCompletion::cost_bound_exhausted);
}

TEST(RaystarCore, UnresolvedInclusiveComparisonFailsWithoutCompletenessClaim) {
  const auto map = makeOpenMap();
  const PlanEndpoint start(5, 5, 5.25, 5.25);
  const PlanEndpoint goal(6, 5, 6.25, 5.25 + std::ldexp(1.0, -40));
  RaystarCore restricted_core(64);

  const auto result =
    restricted_core.plan(map, start, goal, SearchObjective::allWithinCost(1.0), false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.outcome, PlanningOutcome::failed);
  EXPECT_EQ(result.completion, PlanningCompletion::none);
  EXPECT_EQ(result.limit_reached, PlanningLimitReached::none);
  EXPECT_TRUE(result.path_solutions.empty());
  EXPECT_NE(result.message.find("Could not resolve"), std::string::npos);
}

TEST(RaystarCore, CostBoundedPathCapIsAnIncompleteSearchLimit) {
  const auto map = makeSimpleMap();
  RaystarCore core;
  PlanningLimits limits;
  limits.max_cost_bounded_paths = 1;

  const auto limited =
    core.plan(map, 5, 15, 25, 15, SearchObjective::allWithinCost(1000.0), false, limits);
  EXPECT_TRUE(limited.success) << limited.message;
  EXPECT_EQ(limited.outcome, PlanningOutcome::limit_reached);
  EXPECT_EQ(limited.limit_reached, PlanningLimitReached::max_paths);
  EXPECT_EQ(limited.completion, PlanningCompletion::none);
  ASSERT_EQ(limited.path_solutions.size(), 1u);
  EXPECT_NE(limited.message.find("max_cost_bounded_paths=1"), std::string::npos);
}

TEST(RaystarCore, RejectedSearchSupersetPathsDoNotConsumeBoundedPathBudgets) {
  const auto map = makeSimpleMap();
  RaystarCore core;
  const auto baseline = core.plan(map, 5, 14, 25, 15, 2, false);
  ASSERT_TRUE(baseline.success) << baseline.message;
  ASSERT_EQ(baseline.path_solutions.size(), 2u);
  ASSERT_LT(baseline.path_solutions[0].path_cost_, baseline.path_solutions[1].path_cost_);
  const auto& rejected_solution = baseline.path_solutions.front();
  const auto expected_geometry = baseline.path_solutions.back().projectedPath();

  size_t admission_calls = 0;
  const auto admission = [&](const BoundedPathView& path, const StopToken&) {
    ++admission_calls;
    return path.start == rejected_solution.start_ &&
               path.turning_points == rejected_solution.turning_points_ &&
               path.goal == rejected_solution.goal_
             ? BoundedPathAdmission::reject
             : BoundedPathAdmission::accept;
  };
  PlanningLimits limits;
  limits.max_cost_bounded_paths = 1;
  // The retained path fits exactly. The rejected shell must not be charged
  // before the later eligible path is considered.
  limits.max_path_points = expected_geometry.size();
  const auto filtered =
    core.plan(map,
              5,
              14,
              25,
              15,
              SearchObjective::allWithinCost(baseline.path_solutions.back().path_cost_, admission),
              false,
              limits);

  ASSERT_TRUE(filtered.success) << filtered.message;
  EXPECT_EQ(filtered.outcome, PlanningOutcome::complete);
  EXPECT_EQ(filtered.limit_reached, PlanningLimitReached::none);
  EXPECT_TRUE(filtered.completion == PlanningCompletion::cost_bound_exhausted ||
              filtered.completion == PlanningCompletion::frontier_exhausted);
  ASSERT_EQ(filtered.path_solutions.size(), 1u);
  EXPECT_EQ(filtered.path_solutions.front().projectedPath(), expected_geometry);
  EXPECT_GE(admission_calls, 2u);
}

TEST(RaystarCore, BoundedPathAdmissionFailureFailsClosed) {
  const auto map = makeOpenMap();
  RaystarCore core;
  const auto result =
    core.plan(map,
              2,
              2,
              17,
              17,
              SearchObjective::allWithinCost(100.0,
                                             [](const BoundedPathView&, const StopToken&) {
                                               return BoundedPathAdmission::failure_or_stop;
                                             }),
              false);

  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.outcome, PlanningOutcome::failed);
  EXPECT_EQ(result.limit_reached, PlanningLimitReached::none);
  EXPECT_EQ(result.completion, PlanningCompletion::none);
  EXPECT_TRUE(result.path_solutions.empty());
  EXPECT_NE(result.message.find("bounded path admission failed"), std::string::npos);
}

TEST(RaystarCore, BoundedPathAdmissionStopOutranksResourceLimits) {
  const auto map = makeOpenMap();
  RaystarCore core;
  bool cancel_requested = false;
  PlanningLimits limits;
  limits.cancel_requested = [&]() { return cancel_requested; };
  limits.max_cost_bounded_paths = 1;
  limits.max_path_points = 2;
  const auto admission = [&](const BoundedPathView&, const StopToken& stop_token) {
    cancel_requested = true;
    EXPECT_TRUE(stop_token.poll());
    return BoundedPathAdmission::accept;
  };
  const auto result =
    core.plan(map, 2, 2, 17, 17, SearchObjective::allWithinCost(100.0, admission), false, limits);

  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.outcome, PlanningOutcome::limit_reached);
  EXPECT_EQ(result.limit_reached, PlanningLimitReached::cancelled);
  EXPECT_EQ(result.completion, PlanningCompletion::none);
  EXPECT_TRUE(result.path_solutions.empty());
}

TEST(RaystarCore, MultiGoalSingleEntryMatchesSingleGoalBoundedSearch) {
  const auto map = makeSimpleMap();
  const PlanEndpoint start(5, 14, 5.0, 14.0);
  const PlanEndpoint goal(25, 15, 25.0, 15.0);
  constexpr double bound = 40.0;
  RaystarCore core;

  const auto single = core.plan(map, start, goal, SearchObjective::allWithinCost(bound), false);
  const auto multi = core.planToGoalsWithinCosts(map, start, {{goal, bound}}, false);

  ASSERT_EQ(multi.goal_results.size(), 1u);
  const auto& unwrapped = multi.goal_results.front();
  EXPECT_EQ(unwrapped.success, single.success);
  EXPECT_EQ(unwrapped.outcome, single.outcome);
  EXPECT_EQ(unwrapped.limit_reached, single.limit_reached);
  EXPECT_EQ(unwrapped.completion, single.completion);
  ASSERT_EQ(unwrapped.path_solutions.size(), single.path_solutions.size());
  for (size_t i = 0; i < single.path_solutions.size(); ++i) {
    EXPECT_DOUBLE_EQ(unwrapped.path_solutions[i].path_cost_, single.path_solutions[i].path_cost_);
    EXPECT_EQ(unwrapped.path_solutions[i].projectedPath(),
              single.path_solutions[i].projectedPath());
  }
}

TEST(RaystarCore, MultiGoalUsesOneTreeWithIndependentInclusiveBounds) {
  const auto map = makeOpenMap();
  const PlanEndpoint start(2, 2, 2.5, 2.5);
  const PlanEndpoint near_goal(8, 2, 8.5, 2.5);
  const PlanEndpoint far_goal(17, 17, 17.5, 17.5);
  const double near_cost = 6.0;
  const double far_cost = std::hypot(15.0, 15.0);
  RaystarCore core;

  const auto result = core.planToGoalsWithinCosts(
    map, start, {{near_goal, near_cost}, {far_goal, std::nextafter(far_cost, 0.0)}}, false);

  ASSERT_EQ(result.goal_results.size(), 2u);
  EXPECT_EQ(result.outcome, PlanningOutcome::complete);
  ASSERT_EQ(result.goal_results[0].path_solutions.size(), 1u);
  EXPECT_DOUBLE_EQ(result.goal_results[0].path_solutions.front().path_cost_, near_cost);
  EXPECT_TRUE(result.goal_results[0].completion == PlanningCompletion::cost_bound_exhausted ||
              result.goal_results[0].completion == PlanningCompletion::frontier_exhausted);
  EXPECT_FALSE(result.goal_results[1].success);
  EXPECT_EQ(result.goal_results[1].outcome, PlanningOutcome::no_path);
  EXPECT_TRUE(result.goal_results[1].completion == PlanningCompletion::cost_bound_exhausted ||
              result.goal_results[1].completion == PlanningCompletion::frontier_exhausted);
  EXPECT_EQ(result.expanded_nodes, core.getNodes().size());
  EXPECT_EQ(result.expanded_nodes, 1u);
}

TEST(RaystarCore, MultiGoalPreservesOrderAndDuplicateGoalResults) {
  const auto map = makeSimpleMap();
  const PlanEndpoint start(5, 15, 5.5, 15.5);
  const PlanEndpoint upper_goal(25, 14, 25.5, 14.5);
  const PlanEndpoint lower_goal(25, 16, 25.5, 16.5);
  RaystarCore core;

  const auto result = core.planToGoalsWithinCosts(
    map, start, {{upper_goal, 40.0}, {lower_goal, 40.0}, {upper_goal, 40.0}}, false);

  ASSERT_EQ(result.goal_results.size(), 3u);
  EXPECT_EQ(result.goal_results[0].endpoint.position_, upper_goal.position_);
  EXPECT_EQ(result.goal_results[1].endpoint.position_, lower_goal.position_);
  EXPECT_EQ(result.goal_results[2].endpoint.position_, upper_goal.position_);
  ASSERT_EQ(result.goal_results[0].path_solutions.size(),
            result.goal_results[2].path_solutions.size());
  for (size_t i = 0; i < result.goal_results[0].path_solutions.size(); ++i) {
    EXPECT_DOUBLE_EQ(result.goal_results[0].path_solutions[i].path_cost_,
                     result.goal_results[2].path_solutions[i].path_cost_);
    EXPECT_EQ(result.goal_results[0].path_solutions[i].projectedPath(),
              result.goal_results[2].path_solutions[i].projectedPath());
  }
}

TEST(RaystarCore, MultiGoalPathSetsMatchIndependentSearchesAndGoalPermutation) {
  const auto map = makeSimpleMap();
  const PlanEndpoint start(5, 15, 5.5, 15.5);
  const PlanEndpoint upper_goal(25, 14, 25.5, 14.5);
  const PlanEndpoint lower_goal(25, 16, 25.5, 16.5);
  constexpr double bound = 45.0;
  RaystarCore core;

  const auto shared =
    core.planToGoalsWithinCosts(map, start, {{upper_goal, bound}, {lower_goal, bound}}, false);
  ASSERT_EQ(shared.goal_results.size(), 2u);
  const auto upper_independent =
    core.plan(map, start, upper_goal, SearchObjective::allWithinCost(bound), false);
  const auto lower_independent =
    core.plan(map, start, lower_goal, SearchObjective::allWithinCost(bound), false);

  const auto expect_same_paths = [](const std::vector<PathSolution>& lhs,
                                    const std::vector<PathSolution>& rhs) {
    ASSERT_EQ(lhs.size(), rhs.size());
    for (size_t i = 0; i < lhs.size(); ++i) {
      EXPECT_DOUBLE_EQ(lhs[i].path_cost_, rhs[i].path_cost_);
      EXPECT_EQ(lhs[i].projectedPath(), rhs[i].projectedPath());
    }
  };
  expect_same_paths(shared.goal_results[0].path_solutions, upper_independent.path_solutions);
  expect_same_paths(shared.goal_results[1].path_solutions, lower_independent.path_solutions);

  const auto permuted =
    core.planToGoalsWithinCosts(map, start, {{lower_goal, bound}, {upper_goal, bound}}, false);
  ASSERT_EQ(permuted.goal_results.size(), 2u);
  expect_same_paths(shared.goal_results[0].path_solutions, permuted.goal_results[1].path_solutions);
  expect_same_paths(shared.goal_results[1].path_solutions, permuted.goal_results[0].path_solutions);
}

TEST(RaystarCore, MultiGoalSeparatesReachableAndUnreachableGoals) {
  auto map = makeOpenMap();
  for (unsigned int y = 0; y < map.height; ++y) map.data[y * map.width + 10] = 1;
  const PlanEndpoint start(2, 10, 2.5, 10.5);
  const PlanEndpoint reachable_goal(8, 10, 8.5, 10.5);
  const PlanEndpoint unreachable_goal(17, 10, 17.5, 10.5);
  RaystarCore core;

  const auto result = core.planToGoalsWithinCosts(
    map, start, {{reachable_goal, 20.0}, {unreachable_goal, 20.0}}, false);

  ASSERT_EQ(result.goal_results.size(), 2u);
  EXPECT_EQ(result.outcome, PlanningOutcome::complete);
  EXPECT_TRUE(result.goal_results[0].success) << result.goal_results[0].message;
  ASSERT_EQ(result.goal_results[0].path_solutions.size(), 1u);
  EXPECT_FALSE(result.goal_results[1].success);
  EXPECT_EQ(result.goal_results[1].outcome, PlanningOutcome::no_path);
  EXPECT_EQ(result.goal_results[1].completion, PlanningCompletion::frontier_exhausted);
  EXPECT_NE(result.goal_results[1].message.find("No path exists"), std::string::npos);
}

TEST(RaystarCore, MultiGoalRejectsInvalidGoalCountAndBound) {
  const auto map = makeOpenMap();
  const PlanEndpoint start(2, 2, 2.5, 2.5);
  const PlanEndpoint goal(17, 17, 17.5, 17.5);
  RaystarCore core;

  const auto empty = core.planToGoalsWithinCosts(map, start, {}, false);
  EXPECT_EQ(empty.outcome, PlanningOutcome::invalid_request);
  EXPECT_TRUE(core.getNodes().empty());

  const auto invalid_bound = core.planToGoalsWithinCosts(map, start, {{goal, -1.0}}, false);
  EXPECT_EQ(invalid_bound.outcome, PlanningOutcome::invalid_request);
  ASSERT_EQ(invalid_bound.goal_results.size(), 1u);
  EXPECT_EQ(invalid_bound.goal_results.front().outcome, PlanningOutcome::invalid_request);

  PlanningLimits limits;
  limits.max_multi_goal_count = 1;
  const auto too_many =
    core.planToGoalsWithinCosts(map, start, {{goal, 30.0}, {goal, 30.0}}, false, limits);
  EXPECT_EQ(too_many.outcome, PlanningOutcome::invalid_request);
  EXPECT_NE(too_many.message.find("max_multi_goal_count=1"), std::string::npos);
}

TEST(RaystarCore, MultiGoalZeroBoundReturnsOnlyTheIdentityPath) {
  const auto map = makeOpenMap();
  const PlanEndpoint start(2, 2, 2.5, 2.5);
  const PlanEndpoint other_goal(3, 2, 3.5, 2.5);
  RaystarCore core;

  const auto result =
    core.planToGoalsWithinCosts(map, start, {{start, 0.0}, {other_goal, 0.0}}, false);

  ASSERT_EQ(result.goal_results.size(), 2u);
  EXPECT_EQ(result.outcome, PlanningOutcome::complete);

  const auto& identity = result.goal_results[0];
  EXPECT_TRUE(identity.success) << identity.message;
  EXPECT_EQ(identity.outcome, PlanningOutcome::complete);
  ASSERT_EQ(identity.path_solutions.size(), 1u);
  EXPECT_DOUBLE_EQ(identity.path_solutions.front().path_cost_, 0.0);
  EXPECT_EQ(identity.path_solutions.front().projectedPath(),
            std::vector<Point2d>({start.position_, start.position_}));

  const auto& unreachable_within_zero = result.goal_results[1];
  EXPECT_FALSE(unreachable_within_zero.success);
  EXPECT_EQ(unreachable_within_zero.outcome, PlanningOutcome::no_path);
  EXPECT_TRUE(unreachable_within_zero.path_solutions.empty());
}

TEST(RaystarCore, MultiGoalPathCapIsPerGoalAndDoesNotStopOtherGoals) {
  const auto map = makeSimpleMap();
  const PlanEndpoint start(5, 15, 5.5, 15.5);
  const PlanEndpoint obstacle_goal(25, 15, 25.5, 15.5);
  const PlanEndpoint same_side_goal(6, 15, 6.5, 15.5);
  PlanningLimits limits;
  limits.max_cost_bounded_paths = 1;
  RaystarCore core;

  const auto result = core.planToGoalsWithinCosts(
    map, start, {{obstacle_goal, 1000.0}, {same_side_goal, 2.0}}, false, limits);

  ASSERT_EQ(result.goal_results.size(), 2u);
  EXPECT_EQ(result.outcome, PlanningOutcome::limit_reached);
  EXPECT_EQ(result.limit_reached, PlanningLimitReached::max_paths);
  EXPECT_EQ(result.goal_results[0].outcome, PlanningOutcome::limit_reached);
  EXPECT_EQ(result.goal_results[0].limit_reached, PlanningLimitReached::max_paths);
  ASSERT_EQ(result.goal_results[0].path_solutions.size(), 1u);
  EXPECT_EQ(result.goal_results[1].outcome, PlanningOutcome::complete);
  ASSERT_EQ(result.goal_results[1].path_solutions.size(), 1u);
}

TEST(RaystarCore, MultiGoalGlobalPathPointLimitMarksEveryUnfinishedGoal) {
  const auto map = makeOpenMap();
  const PlanEndpoint start(2, 2, 2.5, 2.5);
  const PlanEndpoint first_goal(8, 2, 8.5, 2.5);
  const PlanEndpoint second_goal(2, 8, 2.5, 8.5);
  PlanningLimits limits;
  limits.max_path_points = 3;
  RaystarCore core;

  const auto result = core.planToGoalsWithinCosts(
    map, start, {{first_goal, 10.0}, {second_goal, 10.0}}, false, limits);

  EXPECT_EQ(result.outcome, PlanningOutcome::limit_reached);
  EXPECT_EQ(result.limit_reached, PlanningLimitReached::max_path_points);
  ASSERT_EQ(result.goal_results.size(), 2u);
  for (const auto& goal_result : result.goal_results) {
    EXPECT_EQ(goal_result.outcome, PlanningOutcome::limit_reached);
    EXPECT_EQ(goal_result.limit_reached, PlanningLimitReached::max_path_points);
  }
  EXPECT_EQ(
    result.goal_results[0].path_solutions.size() + result.goal_results[1].path_solutions.size(),
    1u);
}

TEST(RaystarCore, MultiGoalImmediateCancellationIsStructuredAndRecoverable) {
  const auto map = makeOpenMap();
  const PlanEndpoint start(2, 2, 2.5, 2.5);
  const PlanEndpoint goal(17, 17, 17.5, 17.5);
  PlanningLimits limits;
  limits.cancel_requested = []() { return true; };
  RaystarCore core;

  const auto cancelled =
    core.planToGoalsWithinCosts(map, start, {{goal, 30.0}, {goal, 40.0}}, false, limits);
  EXPECT_EQ(cancelled.outcome, PlanningOutcome::limit_reached);
  EXPECT_EQ(cancelled.limit_reached, PlanningLimitReached::cancelled);
  ASSERT_EQ(cancelled.goal_results.size(), 2u);
  for (const auto& goal_result : cancelled.goal_results) {
    EXPECT_EQ(goal_result.outcome, PlanningOutcome::limit_reached);
    EXPECT_EQ(goal_result.limit_reached, PlanningLimitReached::cancelled);
  }
  EXPECT_TRUE(core.getNodes().empty());

  const auto recovered = core.planToGoalsWithinCosts(map, start, {{goal, 30.0}}, false);
  EXPECT_EQ(recovered.outcome, PlanningOutcome::complete);
  ASSERT_EQ(recovered.goal_results.size(), 1u);
  EXPECT_TRUE(recovered.goal_results.front().success);
}

TEST(RaystarCore, RejectsInvalidSearchObjectivesAndBoundedPathLimit) {
  const auto map = makeOpenMap();
  RaystarCore core;

  const std::vector<SearchObjective> invalid_objectives = {
    SearchObjective{SearchMode::all_within_cost, 1, 10.0},
    SearchObjective::allWithinCost(-1.0),
    SearchObjective::allWithinCost(std::numeric_limits<double>::infinity()),
    SearchObjective::allWithinCost(std::numeric_limits<double>::quiet_NaN()),
    SearchObjective{SearchMode::top_k, 1, 10.0}};
  for (const auto& objective : invalid_objectives) {
    const auto result = core.plan(map, 2, 2, 17, 17, objective, false);
    EXPECT_FALSE(result.success);
    EXPECT_EQ(result.outcome, PlanningOutcome::invalid_request);
    EXPECT_TRUE(result.path_solutions.empty());
    EXPECT_TRUE(core.getNodes().empty());
  }

  const auto zero_bound = core.plan(map, 2, 2, 2, 2, SearchObjective::allWithinCost(0.0), false);
  ASSERT_TRUE(zero_bound.success);
  ASSERT_EQ(zero_bound.path_solutions.size(), 1u);
  EXPECT_DOUBLE_EQ(zero_bound.path_solutions.front().path_cost_, 0.0);

  PlanningLimits limits;
  limits.max_cost_bounded_paths = 0;
  const auto invalid_limit =
    core.plan(map, 2, 2, 17, 17, SearchObjective::allWithinCost(100.0), false, limits);
  EXPECT_FALSE(invalid_limit.success);
  EXPECT_NE(invalid_limit.message.find("max_cost_bounded_paths"), std::string::npos);
}

TEST(RaystarCore, ContinuousEndpointsRemainDistinctFromTheirCells) {
  static_assert(std::is_same_v<decltype(Node::seed_), std::pair<int, int>>);
  static_assert(std::is_same_v<decltype(Child::c_), std::pair<int, int>>);

  const auto map = makeOpenMap();
  RaystarCore core;
  const Point2d goal{17.75, 16.5};
  const auto first = core.plan(map, Point2d{2.25, 3.25}, goal, 1, false);
  ASSERT_TRUE(first.success) << first.message;
  ASSERT_EQ(first.path_solutions.size(), 1u);
  const auto first_points = first.path_solutions.front().projectedPath();
  ASSERT_GE(first_points.size(), 2u);
  EXPECT_DOUBLE_EQ(first_points.front().first, 2.25);
  EXPECT_DOUBLE_EQ(first_points.front().second, 3.25);
  EXPECT_DOUBLE_EQ(first_points.back().first, goal.first);
  EXPECT_DOUBLE_EQ(first_points.back().second, goal.second);
  EXPECT_NEAR(first.path_solutions.front().path_cost_,
              std::hypot(goal.first - 2.25, goal.second - 3.25),
              1e-9);

  const auto second = core.plan(map, Point2d{2.75, 3.75}, goal, 1, false);
  ASSERT_TRUE(second.success) << second.message;
  ASSERT_EQ(second.path_solutions.size(), 1u);
  const auto second_points = second.path_solutions.front().projectedPath();
  EXPECT_DOUBLE_EQ(second_points.front().first, 2.75);
  EXPECT_DOUBLE_EQ(second_points.front().second, 3.75);
  EXPECT_NE(second.path_solutions.front().path_cost_, first.path_solutions.front().path_cost_);
  for (size_t index = 1; index + 1 < second_points.size(); ++index) {
    EXPECT_DOUBLE_EQ(second_points[index].first, std::round(second_points[index].first));
    EXPECT_DOUBLE_EQ(second_points[index].second, std::round(second_points[index].second));
  }
}

TEST(RaystarCore, ContinuousEndpointValidationReportsNonFiniteCoordinates) {
  RaystarCore core;
  const auto result = core.plan(makeOpenMap(),
                                Point2d{std::numeric_limits<double>::quiet_NaN(), 2.5},
                                Point2d{17.5, 16.5},
                                1,
                                false);
  EXPECT_FALSE(result.success);
  EXPECT_NE(result.message.find("start"), std::string::npos) << result.message;
  EXPECT_NE(result.message.find("finite"), std::string::npos) << result.message;
}

TEST(RaystarCore, StartEqualsGoalReturnsZeroLengthContinuousPath) {
  RaystarCore core;
  const Point2d endpoint{5.25, 5.75};
  const auto result = core.plan(makeOpenMap(), endpoint, endpoint, 1, false);
  ASSERT_TRUE(result.success) << result.message;
  ASSERT_EQ(result.path_solutions.size(), 1u);
  const auto& solution = result.path_solutions.front();
  EXPECT_DOUBLE_EQ(solution.path_cost_, 0.0);
  EXPECT_EQ(solution.start_, endpoint);
  EXPECT_EQ(solution.goal_, endpoint);
  EXPECT_EQ(solution.turning_points_.size(), 0u);
  ASSERT_EQ(solution.projectedPath().size(), 2u);
  EXPECT_EQ(solution.projectedPath().front(), endpoint);
  EXPECT_EQ(solution.projectedPath().back(), endpoint);
}

TEST(PathSolution, KeepsLegacyIntegerProjectionConstructor) {
  const std::vector<std::pair<int, int>> legacy_path{{2, 3}, {4, 5}, {8, 9}};
  const PathSolution solution(legacy_path, 12.0, {0, 2});
  EXPECT_EQ(solution.path_, legacy_path);
  EXPECT_EQ(solution.start_, std::make_pair(2.0, 3.0));
  EXPECT_EQ(solution.goal_, std::make_pair(8.0, 9.0));
  ASSERT_EQ(solution.turning_points_.size(), 1u);
  EXPECT_EQ(solution.turning_points_.front(), std::make_pair(4, 5));
}

TEST(RaystarCore, RejectsContinuousEndpointOnFinalBoundary) {
  RaystarCore core;
  const auto result = core.plan(makeOpenMap(), Point2d{2.5, 1.0}, Point2d{17.5, 10.5}, 1, false);
  EXPECT_FALSE(result.success);
  EXPECT_NE(result.message.find("start"), std::string::npos) << result.message;
  EXPECT_NE(result.message.find("boundary"), std::string::npos) << result.message;
  EXPECT_TRUE(result.path_solutions.empty());
}

TEST(RaystarCore, RejectsOccupiedAndForcedBorderEndpointsSymmetrically) {
  const auto map_with_occupied_cell = [](int x, int y) {
    auto map = makeOpenMap();
    map.data[static_cast<size_t>(y) * map.width + static_cast<size_t>(x)] = 1;
    return map;
  };

  struct InvalidEndpointCase {
    const char* name;
    GridMap map;
    Point2d start;
    Point2d goal;
    const char* expected_message;
  };

  const auto open_map = makeOpenMap();
  ASSERT_EQ(open_map.data[10 * open_map.width], 0u);
  ASSERT_EQ(open_map.data[10 * open_map.width + open_map.width - 1], 0u);
  const std::vector<InvalidEndpointCase> invalid_cases = {
    {"occupied start",
     map_with_occupied_cell(5, 5),
     Point2d{5.5, 5.5},
     Point2d{17.5, 17.5},
     "Invalid start: corresponding occupancy cell is occupied"},
    {"occupied goal",
     map_with_occupied_cell(15, 15),
     Point2d{2.5, 2.5},
     Point2d{15.5, 15.5},
     "Invalid goal: corresponding occupancy cell is occupied"},
    {"forced-border start",
     open_map,
     Point2d{0.5, 10.5},
     Point2d{17.5, 10.5},
     "Invalid start: corresponding cell is occupied or on the map boundary"},
    {"forced-border goal",
     open_map,
     Point2d{2.5, 10.5},
     Point2d{19.5, 10.5},
     "Invalid goal: corresponding cell is occupied or on the map boundary"},
  };

  RaystarCore core;
  for (const auto& invalid_case : invalid_cases) {
    SCOPED_TRACE(invalid_case.name);
    const auto result =
      core.plan(invalid_case.map, invalid_case.start, invalid_case.goal, 1, false);
    EXPECT_FALSE(result.success);
    EXPECT_EQ(result.outcome, PlanningOutcome::invalid_request);
    EXPECT_EQ(result.limit_reached, PlanningLimitReached::none);
    EXPECT_EQ(result.message, invalid_case.expected_message);
    EXPECT_TRUE(result.path_solutions.empty());
    EXPECT_EQ(result.polymap, nullptr);
    EXPECT_TRUE(core.getNodes().empty());

    const auto recovered =
      core.plan(makeOpenMap(), Point2d{2.5, 2.5}, Point2d{17.5, 17.5}, 1, false);
    ASSERT_TRUE(recovered.success) << recovered.message;
    ASSERT_EQ(recovered.outcome, PlanningOutcome::complete);
    ASSERT_EQ(recovered.path_solutions.size(), 1u);
    ASSERT_FALSE(core.getNodes().empty());
  }
}

TEST(RaystarCore, ContinuousGoalAffectsObstacleRouteAndCost) {
  const auto map = makeSimpleMap();
  RaystarCore core;
  const Point2d start{5.25, 15.25};
  const auto first = core.plan(map, start, Point2d{25.25, 15.25}, 1, false);
  ASSERT_TRUE(first.success) << first.message;
  const auto second = core.plan(map, start, Point2d{25.75, 15.75}, 1, false);
  ASSERT_TRUE(second.success) << second.message;
  ASSERT_EQ(first.path_solutions.size(), 1u);
  ASSERT_EQ(second.path_solutions.size(), 1u);
  EXPECT_DOUBLE_EQ(first.path_solutions.front().goal_.first, 25.25);
  EXPECT_DOUBLE_EQ(second.path_solutions.front().goal_.first, 25.75);
  EXPECT_DOUBLE_EQ(second.path_solutions.front().goal_.second, 15.75);
  EXPECT_NE(first.path_solutions.front().path_cost_, second.path_solutions.front().path_cost_);
  const auto path_length = [](const PathSolution& solution) {
    const auto points = solution.projectedPath();
    double length = 0.0;
    for (size_t i = 1; i < points.size(); ++i)
      length +=
        std::hypot(points[i].first - points[i - 1].first, points[i].second - points[i - 1].second);
    return length;
  };
  EXPECT_NEAR(
    first.path_solutions.front().path_cost_, path_length(first.path_solutions.front()), 1e-9);
  EXPECT_NEAR(
    second.path_solutions.front().path_cost_, path_length(second.path_solutions.front()), 1e-9);
  for (const auto& point : second.path_solutions.front().turning_points_) {
    EXPECT_EQ(point.first, static_cast<int>(point.first));
    EXPECT_EQ(point.second, static_cast<int>(point.second));
  }
}

TEST(RaystarCore, PlanAroundObstacle) {
  auto map = makeSimpleMap();
  RaystarCore core;
  auto result = core.plan(map, 5, 15, 25, 15, 1, false);

  EXPECT_TRUE(result.success);
  EXPECT_GE(result.path_solutions.size(), 1u);
}

TEST(RaystarCore, UpsRemovesCommonTetherPrefixAndCrossesTriangleInteriors) {
  RaystarCore core;
  const Point2d base{5.5, 15.5};
  const auto plan = core.plan(makeSimpleMap(), base, Point2d{25.5, 15.5}, 1, false);
  ASSERT_TRUE(plan.success) << plan.message;
  ASSERT_NE(plan.polymap, nullptr);
  const PathSolution first(
    base, std::vector<std::pair<int, int>>{{10, 10}, {20, 10}}, Point2d{25.5, 15.5}, 0.0, {});
  const PathSolution second(
    base, std::vector<std::pair<int, int>>{{10, 10}, {20, 10}}, Point2d{24.5, 13.5}, 0.0, {});

  const auto transition = RaystarCore::shortenWithinHomotopy(*plan.polymap, first, second);

  ASSERT_TRUE(transition) << transition.message;
  ASSERT_EQ(transition.path.size(), 2u);
  EXPECT_EQ(transition.path.front(), first.goal_);
  EXPECT_EQ(transition.path.back(), second.goal_);
  EXPECT_NEAR(transition.path_cost, std::hypot(1.0, 2.0), 1.0e-12);
}

TEST(RaystarCore, UpsPreservesAConfigurationChangingLoopAtOneRobotPose) {
  RaystarCore core;
  const auto plan = core.plan(makeSimpleMap(), Point2d{5.5, 15.5}, Point2d{25.5, 15.5}, 2, false);
  ASSERT_TRUE(plan.success) << plan.message;
  ASSERT_EQ(plan.path_solutions.size(), 2u);
  ASSERT_NE(plan.polymap, nullptr);

  const auto transition = RaystarCore::shortenWithinHomotopy(
    *plan.polymap, plan.path_solutions[0], plan.path_solutions[1]);

  ASSERT_TRUE(transition) << transition.message;
  ASSERT_GE(transition.path.size(), 2u);
  EXPECT_EQ(transition.path.front(), plan.path_solutions[0].goal_);
  EXPECT_EQ(transition.path.back(), plan.path_solutions[1].goal_);
  EXPECT_EQ(transition.path.front(), transition.path.back());
  EXPECT_GT(transition.path_cost, 0.0);
  EXPECT_TRUE(transition.homotopy_preserved);
  std::string witness_error;
  EXPECT_TRUE(sameReducedDirectedPortalWitness(
    transition.corridor, transition.output_corridor, &witness_error))
    << witness_error;
  EXPECT_EQ(transition.corridor.triangle_occurrences,
            transition.output_corridor.triangle_occurrences);
  bool repeats_a_triangle = false;
  for (size_t first = 0; first < transition.corridor.triangle_occurrences.size(); ++first) {
    for (size_t second = first + 1; second < transition.corridor.triangle_occurrences.size();
         ++second) {
      repeats_a_triangle = repeats_a_triangle || transition.corridor.triangle_occurrences[first] ==
                                                   transition.corridor.triangle_occurrences[second];
    }
  }
  EXPECT_TRUE(repeats_a_triangle)
    << "A lifted winding sleeve must retain repeated geometric triangle occurrences";
}

TEST(RaystarCore, UpsBatchEvaluatesDirectedPairsAndIdentityTransitions) {
  RaystarCore core;
  const auto plan = core.plan(makeSimpleMap(), Point2d{5.5, 15.5}, Point2d{25.5, 15.5}, 2, false);
  ASSERT_TRUE(plan.success) << plan.message;
  ASSERT_EQ(plan.path_solutions.size(), 2u);
  ASSERT_NE(plan.polymap, nullptr);
  const std::vector<ConfigurationTransitionPair> pairs{{0, 0}, {0, 1}, {1, 0}};

  const auto batch =
    RaystarCore::shortenConfigurationTransitions(*plan.polymap, plan.path_solutions, pairs);

  ASSERT_TRUE(batch) << batch.message;
  ASSERT_EQ(batch.transitions.size(), pairs.size());
  EXPECT_EQ(batch.transitions[0].pair.from_configuration, 0u);
  ASSERT_TRUE(batch.transitions[0].shortening) << batch.transitions[0].shortening.message;
  EXPECT_DOUBLE_EQ(batch.transitions[0].shortening.path_cost, 0.0);
  ASSERT_TRUE(batch.transitions[1].shortening) << batch.transitions[1].shortening.message;
  ASSERT_TRUE(batch.transitions[2].shortening) << batch.transitions[2].shortening.message;
  EXPECT_NEAR(
    batch.transitions[1].shortening.path_cost, batch.transitions[2].shortening.path_cost, 1.0e-12);
}

TEST(RaystarCore, UpsBatchRejectsPairIndicesAtomicallyAndSupportsCancellation) {
  RaystarCore core;
  const auto plan = core.plan(makeOpenMap(), Point2d{2.5, 2.5}, Point2d{17.5, 17.5}, 1, false);
  ASSERT_TRUE(plan.success) << plan.message;
  ASSERT_NE(plan.polymap, nullptr);

  const auto invalid = RaystarCore::shortenConfigurationTransitions(
    *plan.polymap, plan.path_solutions, {{0, 0}, {0, 1}});
  EXPECT_EQ(invalid.status, TransitionBatchStatus::invalid_request);
  EXPECT_TRUE(invalid.transitions.empty());

  const StopToken stop([]() { return true; });
  const auto canceled = RaystarCore::shortenConfigurationTransitions(
    *plan.polymap, plan.path_solutions, {{0, 0}}, stop);
  EXPECT_EQ(canceled.status, TransitionBatchStatus::stopped);
  EXPECT_TRUE(canceled.transitions.empty());
}

TEST(RaystarCore, UpsRejectsConfigurationsWithDifferentTetherBases) {
  RaystarCore core;
  const auto plan = core.plan(makeOpenMap(), Point2d{2.5, 2.5}, Point2d{17.5, 17.5}, 1, false);
  ASSERT_TRUE(plan.success) << plan.message;
  ASSERT_NE(plan.polymap, nullptr);
  const PathSolution first(Point2d{2.5, 2.5}, {}, Point2d{4.5, 4.5}, 0.0, {});
  const PathSolution second(Point2d{3.5, 2.5}, {}, Point2d{5.5, 4.5}, 0.0, {});

  const auto transition = RaystarCore::shortenWithinHomotopy(*plan.polymap, first, second);

  EXPECT_EQ(transition.status, HomotopyShorteningStatus::invalid_reference);
  EXPECT_NE(transition.message.find("base"), std::string::npos);
}

TEST(RaystarCore, PathSolutionsRetainCompleteDeepAncestorChains) {
  const auto map = makeSimpleMap();
  RaystarCore core;
  const auto result = core.plan(map, 5, 15, 25, 15, 2, false);

  ASSERT_TRUE(result.success) << result.message;
  ASSERT_EQ(result.outcome, PlanningOutcome::complete);
  ASSERT_EQ(result.path_solutions.size(), 2u);
  const auto& nodes = core.getNodes();
  ASSERT_FALSE(nodes.empty());

  for (size_t solution_index = 0; solution_index < result.path_solutions.size(); ++solution_index) {
    SCOPED_TRACE("solution " + std::to_string(solution_index));
    const auto& solution = result.path_solutions[solution_index];
    const auto& chain = solution.path_node_index_;
    ASSERT_GT(chain.size(), 2u);
    ASSERT_EQ(chain.size(), solution.turning_points_.size() + 1u);
    EXPECT_EQ(chain.front(), 0);

    std::vector<bool> seen(nodes.size(), false);
    for (size_t chain_position = 0; chain_position < chain.size(); ++chain_position) {
      const int node_index = chain[chain_position];
      ASSERT_GE(node_index, 0);
      ASSERT_LT(static_cast<size_t>(node_index), nodes.size());
      EXPECT_FALSE(seen[static_cast<size_t>(node_index)])
        << "ancestor index repeated at chain position " << chain_position;
      seen[static_cast<size_t>(node_index)] = true;

      const auto& node = nodes[static_cast<size_t>(node_index)];
      EXPECT_EQ(node.index(), node_index);
      if (chain_position == 0) {
        EXPECT_EQ(node.parentIndex(), -1);
      } else {
        EXPECT_EQ(node.parentIndex(), chain[chain_position - 1]);
      }
    }

    const int leaf_index = chain.back();
    const auto& leaf = nodes[static_cast<size_t>(leaf_index)];
    EXPECT_EQ(leaf.localShortestPath(), solution.turning_points_);

    // Independently reconstruct root-to-leaf ancestry from parent links.  An
    // invalid cycle is bounded by nodes.size() and fails instead of hanging.
    std::vector<int> reconstructed_reverse;
    int cursor = leaf_index;
    while (cursor >= 0 && reconstructed_reverse.size() <= nodes.size()) {
      ASSERT_LT(static_cast<size_t>(cursor), nodes.size());
      reconstructed_reverse.push_back(cursor);
      cursor = nodes[static_cast<size_t>(cursor)].parentIndex();
    }
    ASSERT_LE(reconstructed_reverse.size(), nodes.size());
    std::reverse(reconstructed_reverse.begin(), reconstructed_reverse.end());
    EXPECT_EQ(reconstructed_reverse, chain);
  }
}

TEST(RaystarCore, EveryExpandedNodeRetainsItsCompleteAncestorPrefix) {
  RaystarCore core;
  const auto result =
    core.plan(makeSimpleMap(), Point2d{5.25, 15.25}, Point2d{25.75, 15.75}, 2, false);

  ASSERT_TRUE(result.success) << result.message;
  const auto& nodes = core.getNodes();
  ASSERT_GT(nodes.size(), 2u);

  for (size_t node_index = 0; node_index < nodes.size(); ++node_index) {
    SCOPED_TRACE("node " + std::to_string(node_index));
    const auto& node = nodes[node_index];
    const auto& prefix = node.pathNodeIndices();
    ASSERT_FALSE(prefix.empty());
    ASSERT_EQ(prefix.back(), static_cast<int>(node_index));
    ASSERT_EQ(prefix.front(), 0);
    ASSERT_LE(prefix.size(), nodes.size());

    std::vector<bool> seen(nodes.size(), false);
    for (size_t depth = 0; depth < prefix.size(); ++depth) {
      const int ancestor_index = prefix[depth];
      ASSERT_GE(ancestor_index, 0);
      ASSERT_LT(static_cast<size_t>(ancestor_index), nodes.size());
      ASSERT_FALSE(seen[static_cast<size_t>(ancestor_index)]);
      seen[static_cast<size_t>(ancestor_index)] = true;
      const auto& ancestor = nodes[static_cast<size_t>(ancestor_index)];
      if (depth == 0)
        EXPECT_EQ(ancestor.parentIndex(), -1);
      else
        EXPECT_EQ(ancestor.parentIndex(), prefix[depth - 1]);
    }
  }
}

TEST(RaystarCore, NoPathToGoal) {
  auto map = makeOpenMap();

  for (unsigned int y = 0; y < 20; ++y) map.data[y * 20 + 10] = 1;

  RaystarCore core;
  auto result = core.plan(map, 5, 10, 15, 10, 1, false);

  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.outcome, PlanningOutcome::no_path);
  EXPECT_EQ(result.limit_reached, PlanningLimitReached::none);
}

TEST(RaystarCore, ReportsUnsupportedSharedVertexTopologyBeforePlanning) {
  GridMap map;
  map.width = 7;
  map.height = 7;
  map.resolution = 1.0f;
  map.data.assign(49, 0);
  map.data[2 * map.width + 2] = 1;
  map.data[3 * map.width + 3] = 1;

  RaystarCore core;
  const auto result = core.plan(map, 4, 2, 5, 5, 1, false);

  EXPECT_FALSE(result.success);
  // A failed controlled construction must not expose the partially built
  // geometry object to callers.
  EXPECT_EQ(result.polymap, nullptr);
  EXPECT_NE(result.message.find("Map geometry construction failed"), std::string::npos)
    << result.message;
  EXPECT_NE(result.message.find("Unsupported obstacle topology"), std::string::npos)
    << result.message;
  EXPECT_NE(result.message.find("(3, 3)"), std::string::npos) << result.message;
  EXPECT_EQ(result.limit_reached, PlanningLimitReached::none);
  EXPECT_TRUE(result.path_solutions.empty());
  EXPECT_TRUE(core.getNodes().empty());
}

TEST(RaystarCore, RejectsMalformedGridMap) {
  RaystarCore core;

  auto map = makeOpenMap();
  map.data.pop_back();
  auto result = core.plan(map, 2, 2, 17, 17, 1, false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid map: data size is 399, expected 400");
  EXPECT_TRUE(result.path_solutions.empty());
  EXPECT_EQ(result.polymap, nullptr);
  EXPECT_TRUE(core.getNodes().empty());

  map = makeOpenMap();
  map.data.push_back(0);
  result = core.plan(map, 2, 2, 17, 17, 1, false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid map: data size is 401, expected 400");
  EXPECT_TRUE(core.getNodes().empty());

  map = makeOpenMap();
  map.width = 0;
  map.data.clear();
  result = core.plan(map, 0, 0, 0, 0, 1, false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid map: width and height must be positive");
  EXPECT_TRUE(core.getNodes().empty());

  GridMap too_wide;
  too_wide.width = static_cast<unsigned int>(std::numeric_limits<int>::max()) + 1U;
  too_wide.height = 1;
  result = core.plan(too_wide, 0, 0, 0, 0, 1, false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid map: width and height must fit signed int indexing");
  EXPECT_TRUE(core.getNodes().empty());

  GridMap too_many_cells;
  too_many_cells.width = 46341;
  too_many_cells.height = 46341;
  result = core.plan(too_many_cells, 0, 0, 0, 0, 1, false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid map: cell count must fit signed int indexing");
  EXPECT_TRUE(core.getNodes().empty());
}

TEST(RaystarCore, EnforcesMapCellAdmissionBeforePlanningStateIsAllocated) {
  const auto map = makeOpenMap();
  RaystarCore core;

  const auto baseline = core.plan(map, 2, 2, 17, 17, 1, false);
  ASSERT_TRUE(baseline.success) << baseline.message;
  ASSERT_FALSE(core.getNodes().empty());

  PlanningLimits limits;
  limits.max_map_cells = map.data.size() - 1;
  const auto rejected = core.plan(
    map, PlanEndpoint(2, 2, 2.25, 2.25), PlanEndpoint(17, 17, 17.75, 17.75), 1, false, limits);

  EXPECT_FALSE(rejected.success);
  EXPECT_EQ(rejected.message, "Invalid map: cell count 400 exceeds max_map_cells=399");
  EXPECT_TRUE(rejected.path_solutions.empty());
  EXPECT_EQ(rejected.polymap, nullptr);
  EXPECT_TRUE(core.getNodes().empty());
}

TEST(RaystarCore, MapCellAdmissionAcceptsTheExactConfiguredBoundary) {
  const auto map = makeOpenMap();
  RaystarCore core;
  PlanningLimits limits;
  limits.max_map_cells = map.data.size();

  const auto result = core.plan(map, 2, 2, 17, 17, 1, false, limits);
  EXPECT_TRUE(result.success) << result.message;
  EXPECT_NE(result.polymap, nullptr);
  EXPECT_FALSE(core.getNodes().empty());
}

TEST(RaystarCore, MapByteAdmissionUsesCheckedWorkingSetEstimate) {
  const auto map = makeOpenMap();
  const size_t estimated_bytes = map.data.size() * kEstimatedPlannerMapBytesPerCell;
  RaystarCore core;

  PlanningLimits exact_limits;
  exact_limits.max_map_bytes = estimated_bytes;
  const auto accepted = core.plan(map, 2, 2, 17, 17, 1, false, exact_limits);
  ASSERT_TRUE(accepted.success) << accepted.message;

  PlanningLimits tight_limits;
  tight_limits.max_map_bytes = estimated_bytes - 1;
  const auto rejected = core.plan(map, 2, 2, 17, 17, 1, false, tight_limits);
  EXPECT_FALSE(rejected.success);
  EXPECT_NE(rejected.message.find("max_map_bytes"), std::string::npos) << rejected.message;
  EXPECT_EQ(rejected.polymap, nullptr);
  EXPECT_TRUE(core.getNodes().empty());

  const auto recovered = core.plan(map, 2, 2, 17, 17, 1, false);
  EXPECT_TRUE(recovered.success) << recovered.message;
}

TEST(RaystarCore, RejectsOutOfBoundsCoordinates) {
  const auto map = makeOpenMap();
  RaystarCore core;

  auto result = core.plan(map, -1, 2, 17, 17, 1, false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid start: coordinates (-1, 2) are outside map bounds");

  result = core.plan(map, 20, 2, 17, 17, 1, false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid start: coordinates (20, 2) are outside map bounds");

  result = core.plan(map, 2, 2, 17, -1, 1, false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid goal: coordinates (17, -1) are outside map bounds");

  result = core.plan(map, 2, 2, 17, 20, 1, false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid goal: coordinates (17, 20) are outside map bounds");
  EXPECT_TRUE(core.getNodes().empty());
}

TEST(RaystarCore, RejectsNonPositiveK) {
  const auto map = makeOpenMap();
  RaystarCore core;

  auto result = core.plan(map, 2, 2, 17, 17, 0, false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid K: must be greater than zero");
  EXPECT_TRUE(core.getNodes().empty());

  result = core.plan(map, 2, 2, 17, 17, -3, false);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid K: must be greater than zero");
  EXPECT_TRUE(core.getNodes().empty());
}

TEST(RaystarCore, RejectsInvalidPlanningLimitsAndKAboveMaxK) {
  const auto map = makeOpenMap();
  RaystarCore core;
  PlanningLimits limits;

  limits.max_k = 3;
  auto result = core.plan(map, 2, 2, 17, 17, 4, false, limits);
  EXPECT_FALSE(result.success);
  EXPECT_EQ(result.message, "Invalid K: requested 4 exceeds max_k=3");
  EXPECT_TRUE(core.getNodes().empty());

  limits = PlanningLimits{};
  limits.max_nodes = 0;
  result = core.plan(map, 2, 2, 17, 17, 1, false, limits);
  EXPECT_FALSE(result.success);
  EXPECT_NE(result.message.find("max_nodes"), std::string::npos);
  EXPECT_TRUE(core.getNodes().empty());

  limits = PlanningLimits{};
  limits.max_nodes = static_cast<size_t>(std::numeric_limits<int>::max()) + 1u;
  result = core.plan(map, 2, 2, 17, 17, 1, false, limits);
  EXPECT_FALSE(result.success);
  EXPECT_NE(result.message.find("max_nodes"), std::string::npos);
  EXPECT_TRUE(core.getNodes().empty());

  limits = PlanningLimits{};
  limits.max_map_cells = 0;
  result = core.plan(map, 2, 2, 17, 17, 1, false, limits);
  EXPECT_FALSE(result.success);
  EXPECT_NE(result.message.find("max_map_cells"), std::string::npos);
  EXPECT_TRUE(core.getNodes().empty());

  limits = PlanningLimits{};
  limits.max_map_bytes = 0;
  result = core.plan(map, 2, 2, 17, 17, 1, false, limits);
  EXPECT_FALSE(result.success);
  EXPECT_NE(result.message.find("max_map_bytes"), std::string::npos);
  EXPECT_TRUE(core.getNodes().empty());

  limits = PlanningLimits{};
  limits.max_path_points = 1;
  result = core.plan(map, 2, 2, 17, 17, 1, false, limits);
  EXPECT_FALSE(result.success);
  EXPECT_NE(result.message.find("max_path_points"), std::string::npos);
  EXPECT_TRUE(core.getNodes().empty());

  limits = PlanningLimits{};
  limits.planning_timeout = std::chrono::milliseconds(-1);
  result = core.plan(map, 2, 2, 17, 17, 1, false, limits);
  EXPECT_FALSE(result.success);
  EXPECT_NE(result.message.find("planning timeout"), std::string::npos);
  EXPECT_TRUE(core.getNodes().empty());
}

TEST(RaystarCore, MaxNodesCountsRootAndHasNoOffByOne) {
  const auto map = makeSimpleMap();
  RaystarCore core;

  const auto baseline = core.plan(map, 5, 15, 25, 15, 1, false);
  ASSERT_TRUE(baseline.success) << baseline.message;
  const size_t required_nodes = core.getNodes().size();
  ASSERT_GT(required_nodes, 1u);

  PlanningLimits limits;
  limits.max_nodes = required_nodes - 1;
  const auto limited = core.plan(map, 5, 15, 25, 15, 1, false, limits);
  EXPECT_FALSE(limited.success);
  EXPECT_EQ(limited.limit_reached, PlanningLimitReached::max_nodes);
  EXPECT_TRUE(limited.path_solutions.empty());
  EXPECT_EQ(core.getNodes().size(), limits.max_nodes);
  EXPECT_NE(limited.message.find("max_nodes"), std::string::npos);

  limits.max_nodes = required_nodes;
  const auto exact_limit = core.plan(map, 5, 15, 25, 15, 1, false, limits);
  EXPECT_TRUE(exact_limit.success) << exact_limit.message;
  EXPECT_EQ(exact_limit.limit_reached, PlanningLimitReached::none);
  EXPECT_EQ(core.getNodes().size(), required_nodes);
}

TEST(RaystarCore, MaxNodesOneStillAllowsRootSolution) {
  const auto map = makeOpenMap();
  RaystarCore core;
  PlanningLimits limits;
  limits.max_nodes = 1;

  const auto result = core.plan(map, 2, 2, 17, 17, 1, false, limits);
  EXPECT_TRUE(result.success) << result.message;
  EXPECT_EQ(result.limit_reached, PlanningLimitReached::none);
  EXPECT_EQ(result.path_solutions.size(), 1u);
  EXPECT_EQ(core.getNodes().size(), 1u);
}

TEST(RaystarCore, ImmediatePlanningTimeoutStopsAndNextPlanRecovers) {
  const auto map = makeOpenMap();
  RaystarCore core;
  PlanningLimits limits;
  limits.planning_timeout = std::chrono::milliseconds(0);

  const auto timed_out = core.plan(map, 2, 2, 17, 17, 1, false, limits);
  EXPECT_FALSE(timed_out.success);
  EXPECT_EQ(timed_out.outcome, PlanningOutcome::limit_reached);
  EXPECT_EQ(timed_out.limit_reached, PlanningLimitReached::timeout);
  EXPECT_TRUE(timed_out.path_solutions.empty());
  EXPECT_TRUE(core.getNodes().empty());
  EXPECT_NE(timed_out.message.find("planning timeout"), std::string::npos);

  const auto recovered = core.plan(map, 2, 2, 17, 17, 1, false);
  EXPECT_TRUE(recovered.success) << recovered.message;
  EXPECT_EQ(recovered.limit_reached, PlanningLimitReached::none);
  EXPECT_FALSE(core.getNodes().empty());
}

TEST(RaystarCore, DeepCooperativeCancellationDoesNotCommitScopedNodeAndRecovers) {
  const auto map = makeSimpleMap();
  RaystarCore calibration_core;
  PlanningLimits calibration_limits;
  size_t scoped_phase_poll_count = 0;
  calibration_limits.cancel_requested = [&]() {
    if (calibration_core.getNodes().size() == 1)
      ++scoped_phase_poll_count;
    return false;
  };

  const auto baseline = calibration_core.plan(map, 5, 15, 25, 15, 1, false, calibration_limits);
  ASSERT_TRUE(baseline.success) << baseline.message;
  ASSERT_GT(calibration_core.getNodes().size(), 1u);
  ASSERT_GT(scoped_phase_poll_count, 0u);

  RaystarCore core;
  PlanningLimits cancellation_limits;
  size_t seen_scoped_phase_polls = 0;
  cancellation_limits.cancel_requested = [&]() {
    if (core.getNodes().size() != 1)
      return false;
    return ++seen_scoped_phase_polls >= scoped_phase_poll_count;
  };

  const auto cancelled = core.plan(map, 5, 15, 25, 15, 1, false, cancellation_limits);
  EXPECT_FALSE(cancelled.success);
  EXPECT_EQ(cancelled.outcome, PlanningOutcome::limit_reached);
  EXPECT_EQ(cancelled.limit_reached, PlanningLimitReached::cancelled);
  EXPECT_EQ(seen_scoped_phase_polls, scoped_phase_poll_count);
  EXPECT_TRUE(cancelled.path_solutions.empty());
  EXPECT_NE(cancelled.message.find("cancelled"), std::string::npos);
  ASSERT_EQ(core.getNodes().size(), 1u);
  EXPECT_TRUE(core.getNodes().front().visibility_region_valid_);
  ASSERT_NE(cancelled.polymap, nullptr);
  EXPECT_FALSE(cancelled.polymap->constructionStopped());
  EXPECT_TRUE(cancelled.polymap->solution_exist_);
  EXPECT_TRUE(cancelled.polymap->isCDTReady());
  EXPECT_TRUE(cancelled.polymap->constructionError().empty());

  const auto recovered = core.plan(map, 5, 15, 25, 15, 1, false);
  EXPECT_TRUE(recovered.success) << recovered.message;
  EXPECT_EQ(recovered.limit_reached, PlanningLimitReached::none);
  EXPECT_FALSE(core.getNodes().empty());
}

TEST(RaystarCore, ExternalAtomicCancellationFromAnotherThreadStopsAndRecovers) {
  const auto map = makeSimpleMap();
  RaystarCore core;
  PlanningLimits cancellation_limits;
  auto cancel_requested = std::make_shared<std::atomic_bool>(false);

  std::mutex poll_mutex;
  std::condition_variable poll_condition;
  bool planner_polled = false;
  cancellation_limits.cancel_requested = [&, cancel_requested]() {
    std::unique_lock<std::mutex> lock(poll_mutex);
    planner_polled = true;
    poll_condition.notify_one();
    poll_condition.wait(lock, [&]() { return cancel_requested->load(std::memory_order_relaxed); });
    return true;
  };

  PlanResult cancelled;
  std::thread planning_thread(
    [&]() { cancelled = core.plan(map, 5, 15, 25, 15, 1, false, cancellation_limits); });

  bool observed_poll = false;
  {
    std::unique_lock<std::mutex> lock(poll_mutex);
    observed_poll =
      poll_condition.wait_for(lock, std::chrono::seconds(5), [&]() { return planner_polled; });
  }
  cancel_requested->store(true, std::memory_order_relaxed);
  poll_condition.notify_all();
  planning_thread.join();

  ASSERT_TRUE(observed_poll);
  EXPECT_FALSE(cancelled.success);
  EXPECT_EQ(cancelled.outcome, PlanningOutcome::limit_reached);
  EXPECT_EQ(cancelled.limit_reached, PlanningLimitReached::cancelled);
  EXPECT_TRUE(cancelled.path_solutions.empty());
  EXPECT_NE(cancelled.message.find("cancelled"), std::string::npos);

  const auto recovered = core.plan(map, 5, 15, 25, 15, 1, false);
  EXPECT_TRUE(recovered.success) << recovered.message;
  EXPECT_EQ(recovered.limit_reached, PlanningLimitReached::none);
  EXPECT_FALSE(core.getNodes().empty());
}

TEST(RaystarCore, NodeLimitReturnsFoundPathsAndMarksSearchIncomplete) {
  const auto map = makeSimpleMap();
  RaystarCore core;
  const auto baseline = core.plan(map, 5, 15, 25, 15, 2, false);
  ASSERT_GE(baseline.path_solutions.size(), 2u) << baseline.message;
  ASSERT_FALSE(baseline.path_solutions.front().path_node_index_.empty());

  PlanningLimits limits;
  limits.max_nodes =
    static_cast<size_t>(baseline.path_solutions.front().path_node_index_.back() + 1);
  const auto limited = core.plan(map, 5, 15, 25, 15, 2, false, limits);

  EXPECT_TRUE(limited.success) << limited.message;
  EXPECT_EQ(limited.outcome, PlanningOutcome::limit_reached);
  EXPECT_EQ(limited.limit_reached, PlanningLimitReached::max_nodes);
  EXPECT_FALSE(limited.path_solutions.empty());
  EXPECT_LT(limited.path_solutions.size(), 2u);
  EXPECT_NE(limited.message.find("search is incomplete"), std::string::npos);
  EXPECT_LE(core.getNodes().size(), limits.max_nodes);
  EXPECT_EQ(limited.expanded_nodes, core.getNodes().size());
}

TEST(Node, GenerateChildRejectsDegenerateVisibility) {
  auto map = makeSimpleMap();
  Polymap poly(map, 5, 15, 25, 15);
  ASSERT_TRUE(poly.solution_exist_);

  const std::vector<std::pair<double, double>> empty_visibility;
  const std::vector<std::pair<int, int>> empty_topology;
  Node empty_node(&poly, 0, 5, 5, 0.0, 0.0, empty_visibility, empty_topology);
  empty_node.C_.emplace_back(Child(0, 0, 6, 5, false));

  std::string error;
  EXPECT_FALSE(empty_node.generateChild(&poly, &error));
  EXPECT_FALSE(error.empty());
  EXPECT_TRUE(empty_node.C_.empty());

  const std::vector<std::pair<double, double>> visibility = {{6.0, 5.0}, {7.0, 5.0}};
  const std::vector<std::pair<int, int>> mismatched_topology = {{-1, -1}};
  Node mismatched_node(&poly, 1, 5, 5, 0.0, 0.0, visibility, mismatched_topology);

  error.clear();
  EXPECT_FALSE(mismatched_node.generateChild(&poly, &error));
  EXPECT_FALSE(error.empty());
  EXPECT_TRUE(mismatched_node.C_.empty());
}

TEST(Node, CooperativeStopDoesNotCommitFullVisibilityOrChildren) {
  auto map = makeSimpleMap();
  for (unsigned int x = 0; x < map.width; ++x) {
    map.data[x] = 1;
    map.data[(map.height - 1) * map.width + x] = 1;
  }
  for (unsigned int y = 0; y < map.height; ++y) {
    map.data[y * map.width] = 1;
    map.data[y * map.width + map.width - 1] = 1;
  }
  Polymap poly(map, 5, 15, 25, 15);
  ASSERT_TRUE(poly.isCDTReady()) << poly.constructionError();

  VisibilityRegion visibility;
  std::string error;
  ASSERT_TRUE(poly.getVisibilityRegion(5, 15, visibility, &error)) << error;
  ASSERT_FALSE(visibility.empty());

  Node node(&poly, 0, 5, 15, 0.0, 0.0, visibility);
  ASSERT_TRUE(node.visibility_region_valid_) << node.visibility_region_error_;
  node.setFullVisibilityRegion(visibility);
  const auto original_full_visibility = node.full_visibility_region_;

  VisibilityRegion replacement;
  replacement.reserve(128);
  for (size_t repeat = 0; repeat < 128; ++repeat)
    replacement.emplace_back(visibility[repeat % visibility.size()]);

  size_t projection_full_poll_count = 0;
  const StopToken projection_count_only([&projection_full_poll_count]() {
    ++projection_full_poll_count;
    return false;
  });
  ASSERT_EQ(node.setFullVisibilityRegion(replacement, projection_count_only),
            OperationStatus::success);
  ASSERT_GT(projection_full_poll_count, 0u);
  node.setFullVisibilityRegion(original_full_visibility);

  size_t projection_stopped_poll_count = 0;
  const StopToken projection_stop_at_last(
    [&projection_stopped_poll_count, projection_full_poll_count]() {
      return ++projection_stopped_poll_count >= projection_full_poll_count;
    });
  EXPECT_EQ(node.setFullVisibilityRegion(replacement, projection_stop_at_last),
            OperationStatus::stopped);
  EXPECT_EQ(projection_stopped_poll_count, projection_full_poll_count);
  ASSERT_EQ(node.full_visibility_region_.size(), original_full_visibility.size());
  for (size_t index = 0; index < original_full_visibility.size(); ++index) {
    EXPECT_EQ(exactPoint(node.full_visibility_region_[index]),
              exactPoint(original_full_visibility[index]));
    EXPECT_EQ(node.full_visibility_region_[index].support, original_full_visibility[index].support);
  }

  size_t child_full_poll_count = 0;
  const StopToken child_count_only([&child_full_poll_count]() {
    ++child_full_poll_count;
    return false;
  });
  ASSERT_EQ(node.generateChild(&poly, child_count_only, &error), OperationStatus::success) << error;
  ASSERT_GT(child_full_poll_count, 0u);
  ASSERT_FALSE(node.C_.empty());
  const auto expected_children = node.C_;

  error = "stale error";
  size_t child_stopped_poll_count = 0;
  const StopToken child_stop_at_last([&child_stopped_poll_count, child_full_poll_count]() {
    return ++child_stopped_poll_count >= child_full_poll_count;
  });
  EXPECT_EQ(node.generateChild(&poly, child_stop_at_last, &error), OperationStatus::stopped);
  EXPECT_EQ(child_stopped_poll_count, child_full_poll_count);
  EXPECT_TRUE(node.C_.empty());
  EXPECT_TRUE(error.empty());

  ASSERT_TRUE(node.generateChild(&poly, &error)) << error;
  ASSERT_EQ(node.C_.size(), expected_children.size());
  for (size_t index = 0; index < expected_children.size(); ++index) {
    const auto& actual = node.C_[index];
    const auto& expected = expected_children[index];
    EXPECT_EQ(actual.c_, expected.c_);
    EXPECT_EQ(exactPoint(actual.c_endpoint_), exactPoint(expected.c_endpoint_));
    EXPECT_EQ(actual.c_endpoint_.support, expected.c_endpoint_.support);
    EXPECT_EQ(exactPoint(actual.o_endpoint_), exactPoint(expected.o_endpoint_));
    EXPECT_EQ(actual.o_endpoint_.support, expected.o_endpoint_.support);
    EXPECT_EQ(actual.is_a_left_gap_, expected.is_a_left_gap_);
    EXPECT_DOUBLE_EQ(actual.start_angle_, expected.start_angle_);
    EXPECT_DOUBLE_EQ(actual.end_angle_, expected.end_angle_);
    EXPECT_DOUBLE_EQ(actual.c_gcost_, expected.c_gcost_);
  }
}

TEST(Node, GenerateChildRejectsCoordinatesOutsideIntRange) {
  auto map = makeSimpleMap();
  Polymap poly(map, 5, 15, 25, 15);
  ASSERT_TRUE(poly.solution_exist_);

  const double above_int = static_cast<double>(std::numeric_limits<int>::max()) + 1024.0;
  const double below_int = static_cast<double>(std::numeric_limits<int>::lowest()) - 1024.0;
  const std::vector<std::vector<std::pair<double, double>>> invalid_regions = {
    {{above_int, 5.0}, {above_int + 1024.0, 5.0}}, {{below_int, 5.0}, {below_int - 1024.0, 5.0}}};
  const std::vector<std::pair<int, int>> gap_topology = {{0, 0}, {-1, -1}};

  for (size_t i = 0; i < invalid_regions.size(); ++i) {
    SCOPED_TRACE(i);
    Node node(&poly, static_cast<int>(i), 5, 5, 0.0, 0.0, invalid_regions[i], gap_topology);
    node.C_.emplace_back(Child(static_cast<int>(i), 0, 6, 5, false));

    std::string error;
    EXPECT_FALSE(node.generateChild(&poly, &error));
    EXPECT_NE(error.find("representable int range"), std::string::npos);
    EXPECT_TRUE(node.C_.empty());
  }
}

TEST(Node, LegacyAdapterRejectsNonFiniteCoordinatesBeforeExactConstruction) {
  auto map = makeSimpleMap();
  Polymap poly(map, 5, 15, 25, 15);
  ASSERT_TRUE(poly.solution_exist_);

  const std::vector<double> invalid_coordinates = {std::numeric_limits<double>::quiet_NaN(),
                                                   std::numeric_limits<double>::infinity(),
                                                   -std::numeric_limits<double>::infinity()};
  const std::vector<std::pair<int, int>> topology = {{-1, -1}, {-1, -1}};

  for (size_t i = 0; i < invalid_coordinates.size(); ++i) {
    SCOPED_TRACE(i);
    const std::vector<std::pair<double, double>> visibility = {{invalid_coordinates[i], 5.0},
                                                               {7.0, 5.0}};
    Node node(&poly, static_cast<int>(i), 5, 5, 0.0, 0.0, visibility, topology);
    EXPECT_FALSE(node.visibility_region_valid_);

    std::string error;
    EXPECT_FALSE(node.generateChild(&poly, &error));
    EXPECT_NE(error.find("representable int range"), std::string::npos);
    EXPECT_TRUE(node.C_.empty());
  }
}

TEST(Node, ConstructorsSafelyRejectInvalidSeeds) {
  auto map = makeSimpleMap();
  Polymap poly(map, 5, 15, 25, 15);
  ASSERT_TRUE(poly.solution_exist_);

  const std::vector<std::pair<double, double>> visibility = {{6.0, 5.0}, {7.0, 5.0}};
  const std::vector<std::pair<int, int>> topology = {{-1, -1}, {-1, -1}};

  Node root_node(
    &poly, 0, std::numeric_limits<double>::infinity(), 5.0, 0.0, 0.0, visibility, topology);
  std::string error;
  EXPECT_FALSE(root_node.generateChild(&poly, &error));
  EXPECT_NE(error.find("Node seed"), std::string::npos);
  EXPECT_TRUE(root_node.C_.empty());

  const double below_int = static_cast<double>(std::numeric_limits<int>::lowest()) - 1024.0;
  Node scoped_node(&poly, 1, below_int, 5.0, 0.0, 0.0, 0, visibility, topology);
  error.clear();
  EXPECT_FALSE(scoped_node.generateChild(&poly, &error));
  EXPECT_NE(error.find("Node seed"), std::string::npos);
  EXPECT_TRUE(scoped_node.C_.empty());
}

TEST(Child, ExactEndpointProjectionRoundsCanonicalGridCoordinate) {
  auto map = makeSimpleMap();
  Polymap poly(map, 5, 15, 25, 15);
  ASSERT_TRUE(poly.solution_exist_);
  ASSERT_FALSE(poly.obs_.empty());
  ASSERT_FALSE(poly.obs_[0].ordered_vertices_.empty());

  const auto coordinate = poly.obs_[0].ordered_vertices_[0];
  BoundaryEndpoint near_endpoint;
  near_endpoint.position = {static_cast<double>(coordinate.first),
                            static_cast<double>(coordinate.second)};
  if (coordinate.first > 0)
    near_endpoint.position.first -= 5e-10;
  else
    near_endpoint.position.second -= 5e-10;
  near_endpoint.support = ObstacleVertexId{0, 0};

  std::string error;
  EXPECT_FALSE(poly.validateBoundaryEndpoint(near_endpoint, &error));
  EXPECT_NE(error.find("exact position"), std::string::npos);

  const VisibilityRegion visibility = {near_endpoint, near_endpoint};
  Node node(&poly, 0, 5, 5, 0.0, 0.0, visibility);
  ASSERT_TRUE(node.visibility_region_valid_) << node.visibility_region_error_;
  ASSERT_EQ(node.visibility_region_.size(), 2u);
  EXPECT_EQ(
    node.visibility_region_[0].position,
    (Point2d{static_cast<double>(coordinate.first), static_cast<double>(coordinate.second)}));
  EXPECT_EQ(exactPoint(node.visibility_region_[0]),
            exact_geometry::Point(coordinate.first, coordinate.second));

  Child child(0, 0, node.visibility_region_[0], node.visibility_region_[1], false);
  EXPECT_EQ(child.c_, coordinate);
  EXPECT_TRUE(poly.validateBoundaryEndpoint(child.c_endpoint_, &error)) << error;
}

TEST(Node, GenerateChildIgnoresZeroLengthSamePointDiscontinuity) {
  auto map = makeSimpleMap();
  Polymap poly(map, 5, 15, 25, 15);
  ASSERT_TRUE(poly.solution_exist_);
  ASSERT_FALSE(poly.obs_.empty());
  ASSERT_FALSE(poly.obs_[0].ordered_vertices_.empty());

  const auto coordinate = poly.obs_[0].ordered_vertices_[0];
  const BoundaryEndpoint repeated_endpoint{
    {static_cast<double>(coordinate.first), static_cast<double>(coordinate.second)},
    ObstacleVertexId{0, 0}};
  const VisibilityRegion visibility = {repeated_endpoint, repeated_endpoint};
  Node node(&poly, 0, 5, 5, 0.0, 0.0, visibility);
  ASSERT_TRUE(node.visibility_region_valid_) << node.visibility_region_error_;

  std::string error;
  EXPECT_TRUE(node.generateChild(&poly, &error)) << error;
  EXPECT_TRUE(node.C_.empty());
}

TEST(Node, GenerateChildHandlesUnknownAndRejectsInvalidTopology) {
  auto map = makeSimpleMap();
  Polymap poly(map, 5, 15, 25, 15);
  ASSERT_TRUE(poly.solution_exist_);
  ASSERT_FALSE(poly.obs_.empty());

  const std::vector<std::pair<double, double>> visibility = {{6.0, 5.0}, {7.0, 5.0}};
  const std::vector<std::pair<int, int>> unknown_topology = {{-1, -1}, {-1, -1}};
  Node unknown_node(&poly, 0, 5, 5, 0.0, 0.0, visibility, unknown_topology);

  std::string error;
  EXPECT_TRUE(unknown_node.generateChild(&poly, &error)) << error;
  EXPECT_TRUE(error.empty());
  EXPECT_TRUE(unknown_node.C_.empty());

  const std::vector<std::pair<int, int>> half_invalid_topology = {{-1, 0}, {-1, -1}};
  Node invalid_node(&poly, 1, 5, 5, 0.0, 0.0, visibility, half_invalid_topology);

  error.clear();
  EXPECT_FALSE(invalid_node.generateChild(&poly, &error));
  EXPECT_FALSE(error.empty());
  EXPECT_TRUE(invalid_node.C_.empty());

  const std::vector<std::pair<int, int>> out_of_range_topology = {
    {0, static_cast<int>(poly.obs_[0].ordered_vertices_.size())}, {-1, -1}};
  Node out_of_range_node(&poly, 2, 5, 5, 0.0, 0.0, visibility, out_of_range_topology);

  error.clear();
  EXPECT_FALSE(out_of_range_node.generateChild(&poly, &error));
  EXPECT_FALSE(error.empty());
  EXPECT_TRUE(out_of_range_node.C_.empty());
}

TEST(RaystarCore, ClearsNodesBeforeEveryPlan) {
  auto map = makeOpenMap();
  RaystarCore core;

  auto first_result = core.plan(map, 2, 2, 17, 17, 1, false);
  ASSERT_TRUE(first_result.success);
  ASSERT_FALSE(core.getNodes().empty());

  auto rejected_result = core.plan(map, 2, 2, 17, 17, 0, false);
  EXPECT_FALSE(rejected_result.success);
  EXPECT_TRUE(core.getNodes().empty());

  auto second_result = core.plan(map, 2, 2, 17, 17, 1, false);
  ASSERT_TRUE(second_result.success);
  ASSERT_FALSE(core.getNodes().empty());

  for (unsigned int y = 0; y < map.height; ++y) map.data[y * map.width + 10] = 1;

  auto no_path_result = core.plan(map, 5, 10, 15, 10, 1, false);
  EXPECT_FALSE(no_path_result.success);
  EXPECT_EQ(no_path_result.message, "No path exists between start and goal");
  EXPECT_TRUE(core.getNodes().empty());

  auto recovered_result = core.plan(makeOpenMap(), 2, 2, 17, 17, 1, false);
  EXPECT_TRUE(recovered_result.success);
  EXPECT_FALSE(core.getNodes().empty());

  core.resetSearchState();
  EXPECT_TRUE(core.getNodes().empty());
}

TEST(RaystarCore, MultiplePaths) {
  auto map = makeSimpleMap();
  RaystarCore core;
  auto result = core.plan(map, 5, 15, 25, 15, 5, false);

  EXPECT_TRUE(result.success);
  EXPECT_GE(result.path_solutions.size(), 1u);
}

TEST(RaystarCore, SelfCrossingPruning) {
  auto map = makeSimpleMap();
  RaystarCore core;

  auto result_with = core.plan(map, 5, 15, 25, 15, 3, true);
  auto result_without = core.plan(map, 5, 15, 25, 15, 3, false);

  EXPECT_TRUE(result_with.success);
  EXPECT_TRUE(result_without.success);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
