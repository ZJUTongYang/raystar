#include <gtest/gtest.h>

#include <raystar/polymap.h>

#include "visibility_validation.h"

#include <algorithm>
#include <cctype>
#include <initializer_list>
#include <optional>
#include <string>
#include <utility>

namespace {

using raystar::BoundaryEndpoint;
using raystar::OperationStatus;
using raystar::StopToken;
using raystar::VisibilityRegion;
using raystar::detail::VisibilityBoundaryMode;
using raystar::detail::VisibilityGeometryContext;
using raystar::exact_geometry::Point;

BoundaryEndpoint makeEndpoint(int x, int y) {
  BoundaryEndpoint endpoint;
  endpoint.position = {static_cast<double>(x), static_cast<double>(y)};
  endpoint.exact_position = Point(x, y);
  endpoint.support = std::monostate{};
  return endpoint;
}

VisibilityRegion makeBoundary(std::initializer_list<std::pair<int, int>> vertices) {
  VisibilityRegion boundary;
  boundary.reserve(vertices.size());
  for (const auto& [x, y] : vertices) boundary.emplace_back(makeEndpoint(x, y));
  return boundary;
}

VisibilityGeometryContext closedRootContext(const Point& source) {
  return VisibilityGeometryContext{
    source, VisibilityBoundaryMode::closed_cycle, std::nullopt, std::nullopt};
}

VisibilityGeometryContext openSectorContext(const Point& source,
                                            const Point& start,
                                            const Point& end) {
  return VisibilityGeometryContext{source, VisibilityBoundaryMode::open_sector, start, end};
}

std::string lowerCase(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
    return static_cast<char>(std::tolower(character));
  });
  return value;
}

void expectRejected(const VisibilityRegion& boundary,
                    const VisibilityGeometryContext& context,
                    const std::string& expected_error_fragment) {
  std::string error = "stale error";
  EXPECT_FALSE(raystar::detail::validateVisibilityGeometry(boundary, context, &error));
  ASSERT_FALSE(error.empty());
  EXPECT_NE(lowerCase(error).find(expected_error_fragment), std::string::npos) << error;
}

TEST(VisibilityGeometryValidation, ClosedRootAllowsAdjacentRepeatedEndpoint) {
  const auto boundary = makeBoundary({{-4, -4}, {4, -4}, {4, -4}, {4, 4}, {-4, 4}});
  std::string error = "stale error";

  EXPECT_TRUE(
    raystar::detail::validateVisibilityGeometry(boundary, closedRootContext(Point(0, 0)), &error))
    << error;
  EXPECT_TRUE(error.empty());
}

TEST(VisibilityGeometryValidation, ClosedRootAllowsSameRayBacktrackAndOverlap) {
  // The three vertical points are on the same directed ray from the source.
  // The far -> near -> far pair creates two collinear, oppositely directed
  // boundary segments. This is a legitimate visibility discontinuity, not a
  // proper polygon crossing.
  const auto boundary = makeBoundary({{-4, -4}, {4, -4}, {4, 4}, {0, 6}, {0, 3}, {0, 6}, {-4, 4}});
  std::string error = "stale error";

  EXPECT_TRUE(
    raystar::detail::validateVisibilityGeometry(boundary, closedRootContext(Point(0, 0)), &error))
    << error;
  EXPECT_TRUE(error.empty());
}

TEST(VisibilityGeometryValidation, ClosedRootRejectsSourceEndpoint) {
  const auto boundary = makeBoundary({{-4, -4}, {4, -4}, {0, 0}, {4, 4}, {-4, 4}});

  expectRejected(boundary, closedRootContext(Point(0, 0)), "source");
}

TEST(VisibilityGeometryValidation, ClosedRootRejectsRayOrderRegression) {
  // A three-edge boundary cannot have a proper self-crossing. Its rays are
  // nevertheless emitted as -135, +90, 0 degrees instead of monotonically in
  // the fixed atan2 order.
  const auto boundary = makeBoundary({{-4, -4}, {0, 6}, {4, 0}});

  expectRejected(boundary, closedRootContext(Point(0, 0)), "order");
}

TEST(VisibilityGeometryValidation, ClosedRootRejectsProperCrossing) {
  // Segments 0->1 and 2->3 cross at (0,2), strictly inside both segments.
  // The source is deliberately elsewhere so this exercises the crossing
  // invariant rather than the source-on-boundary invariant.
  const auto boundary = makeBoundary({{-4, 0}, {4, 4}, {-4, 4}, {4, 0}});

  expectRejected(boundary, closedRootContext(Point(0, -10)), "cross");
}

TEST(VisibilityGeometryValidation, OpenSectorAcceptsMatchingAnchors) {
  const Point source(0, 0);
  const Point start(4, 0);
  const Point end(0, 4);
  const auto boundary = makeBoundary({{4, 0}, {4, 4}, {0, 4}});
  std::string error = "stale error";

  EXPECT_TRUE(raystar::detail::validateVisibilityGeometry(
    boundary, openSectorContext(source, start, end), &error))
    << error;
  EXPECT_TRUE(error.empty());
}

TEST(VisibilityGeometryValidation, OpenSectorAllowsSameRayBacktrackAndOverlap) {
  const Point source(10, 10);
  const Point start(6, 6);
  const Point end(6, 14);
  const auto boundary =
    makeBoundary({{6, 6}, {14, 6}, {14, 14}, {10, 16}, {10, 13}, {10, 16}, {6, 14}});
  std::string error = "stale error";

  EXPECT_TRUE(raystar::detail::validateVisibilityGeometry(
    boundary, openSectorContext(source, start, end), &error))
    << error;
  EXPECT_TRUE(error.empty());
}

TEST(VisibilityGeometryValidation, OpenSectorDoesNotInventClosingChord) {
  const Point source(0, 0);
  const Point start(4, 0);
  const Point end(0, 4);
  // The artificial end->start chord crosses the vertical edge at (2, 2), but
  // an open visibility fan has source spokes instead of that chord. Its real
  // conceptual boundary is valid and must therefore be accepted.
  const auto boundary = makeBoundary({{4, 0}, {2, 1}, {2, 4}, {0, 4}});
  std::string error = "stale error";

  EXPECT_TRUE(raystar::detail::validateVisibilityGeometry(
    boundary, openSectorContext(source, start, end), &error))
    << error;
  EXPECT_TRUE(error.empty());
}

TEST(VisibilityGeometryValidation, OpenSectorRejectsSpokeCrossing) {
  const Point source(0, 0);
  const Point start(-1, 1);
  const Point end(1, 1);
  // The wide CCW sector contains every endpoint in order, but start->(2,-1)
  // crosses the nonadjacent end->source spoke strictly inside both edges.
  const auto boundary = makeBoundary({{-1, 1}, {2, -1}, {1, 1}});

  expectRejected(boundary, openSectorContext(source, start, end), "cross");
}

TEST(VisibilityGeometryValidation, OpenSectorRejectsAnchorMismatch) {
  const Point source(0, 0);
  const Point start(4, 0);
  const Point end(0, 4);
  const auto context = openSectorContext(source, start, end);

  const std::vector<VisibilityRegion> mismatched_boundaries = {
    makeBoundary({{2, 0}, {4, 4}, {0, 4}}), makeBoundary({{4, 0}, {4, 4}, {0, 2}})};

  for (size_t index = 0; index < mismatched_boundaries.size(); ++index) {
    SCOPED_TRACE(index);
    expectRejected(mismatched_boundaries[index], context, "anchor");
  }
}

TEST(VisibilityGeometryValidation, OpenSectorRejectsRayOutsideSector) {
  const Point source(0, 0);
  const Point start(4, 0);
  const Point end(0, 4);
  const auto context = openSectorContext(source, start, end);

  const std::vector<VisibilityRegion> out_of_sector_boundaries = {
    makeBoundary({{4, 0}, {-1, 6}, {0, 4}}), makeBoundary({{4, 0}, {-1, -1}, {0, 4}})};

  for (size_t index = 0; index < out_of_sector_boundaries.size(); ++index) {
    SCOPED_TRACE(index);
    expectRejected(out_of_sector_boundaries[index], context, "sector");
  }
}

TEST(VisibilityGeometryValidation, OpenSectorRejectsRayOrderRegression) {
  const Point source(0, 0);
  const Point start(4, 0);
  const Point end(0, 4);
  // All rays lie inside the first-quadrant sector, but the direction regresses
  // from (1, 2) to (4, 1). The conceptual source-spoke closure stays simple,
  // so this reaches the ray-order check rather than the crossing check.
  const auto boundary = makeBoundary({{4, 0}, {1, 2}, {4, 1}, {0, 4}});

  expectRejected(boundary, openSectorContext(source, start, end), "order");
}

TEST(VisibilityGeometryValidation, StopsDuringQuadraticCrossingCheck) {
  constexpr size_t endpoint_count = 256;
  VisibilityRegion boundary;
  boundary.reserve(endpoint_count);
  for (size_t index = 0; index < endpoint_count; ++index) boundary.emplace_back(makeEndpoint(1, 0));

  // Entry plus all pre-crossing linear passes require fewer than 800 polls.
  // Reaching 1100 therefore deterministically stops inside the O(B^2) edge
  // pair scan, while still leaving ample pairs for the predicate to trigger.
  size_t poll_count = 0;
  const StopToken stop_token([&poll_count]() { return ++poll_count >= 1100; });
  std::string error = "stale error";

  EXPECT_EQ(raystar::detail::validateVisibilityGeometry(
              boundary, closedRootContext(Point(0, 0)), stop_token, &error),
            OperationStatus::stopped);
  EXPECT_TRUE(stop_token.requested());
  EXPECT_TRUE(error.empty());

  // A stopped validation has no process-wide effect: the legacy no-stop
  // wrapper retains both its success and failure semantics afterward.
  EXPECT_TRUE(
    raystar::detail::validateVisibilityGeometry(boundary, closedRootContext(Point(0, 0)), &error))
    << error;
  EXPECT_TRUE(error.empty());

  const auto crossing_boundary = makeBoundary({{-4, 0}, {4, 4}, {-4, 4}, {4, 0}});
  const StopToken never_stop;
  EXPECT_EQ(raystar::detail::validateVisibilityGeometry(
              crossing_boundary, closedRootContext(Point(0, -10)), never_stop, &error),
            OperationStatus::failure);
  EXPECT_NE(lowerCase(error).find("cross"), std::string::npos) << error;
}

}  // namespace
