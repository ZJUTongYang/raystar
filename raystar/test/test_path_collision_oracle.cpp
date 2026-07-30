#include "path_collision_oracle.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

using raystar::test_oracle::OracleOptions;
using raystar::test_oracle::SelfIntersectionKind;
using raystar::test_oracle::SelfIntersectionStatus;
using raystar::test_oracle::validateNoSelfIntersection;
using raystar::test_oracle::validatePath;
using raystar::test_oracle::validateSegment;
using raystar::test_oracle::ValidationStatus;

nav_msgs::msg::OccupancyGrid makeGrid(std::uint32_t width = 8,
                                      std::uint32_t height = 8,
                                      float resolution = 1.0F,
                                      double origin_x = 0.0,
                                      double origin_y = 0.0) {
  nav_msgs::msg::OccupancyGrid map;
  map.header.frame_id = "map";
  map.info.width = width;
  map.info.height = height;
  map.info.resolution = resolution;
  map.info.origin.position.x = origin_x;
  map.info.origin.position.y = origin_y;
  map.info.origin.orientation.w = 1.0;
  map.data.assign(static_cast<std::size_t>(width) * height, 0);
  return map;
}

void setCell(nav_msgs::msg::OccupancyGrid& map,
             std::uint32_t x,
             std::uint32_t y,
             std::int8_t value) {
  map.data[static_cast<std::size_t>(y) * map.info.width + x] = value;
}

geometry_msgs::msg::Point gridPoint(const nav_msgs::msg::OccupancyGrid& map, double x, double y) {
  geometry_msgs::msg::Point point;
  point.x = map.info.origin.position.x + x * static_cast<double>(map.info.resolution);
  point.y = map.info.origin.position.y + y * static_cast<double>(map.info.resolution);
  return point;
}

nav_msgs::msg::Path makePath(std::initializer_list<std::pair<double, double>> points) {
  nav_msgs::msg::Path path;
  path.header.frame_id = "map";
  for (const auto& point : points) {
    geometry_msgs::msg::PoseStamped pose;
    pose.header.frame_id = path.header.frame_id;
    pose.pose.position.x = point.first;
    pose.pose.position.y = point.second;
    path.poses.push_back(pose);
  }
  return path;
}

void expectCollision(const nav_msgs::msg::OccupancyGrid& map,
                     const geometry_msgs::msg::Point& start,
                     const geometry_msgs::msg::Point& end,
                     const OracleOptions& options = OracleOptions{}) {
  const auto result = validateSegment(map, start, end, options);
  EXPECT_EQ(result.status, ValidationStatus::kCollision) << result.diagnostic;
  EXPECT_TRUE(result.collisionDetected()) << result.diagnostic;
  EXPECT_FALSE(result.diagnostic.empty());
}

void expectFree(const nav_msgs::msg::OccupancyGrid& map,
                const geometry_msgs::msg::Point& start,
                const geometry_msgs::msg::Point& end,
                const OracleOptions& options = OracleOptions{}) {
  const auto result = validateSegment(map, start, end, options);
  EXPECT_EQ(result.status, ValidationStatus::kCollisionFree) << result.diagnostic;
  EXPECT_TRUE(result.collisionFree()) << result.diagnostic;
}

TEST(PathCollisionOracle, DistinguishesInteriorFromLegalBoundaryAndVertexContact) {
  auto map = makeGrid();
  setCell(map, 3, 3, 100);

  // Ordinary and diagonal crossings enter the occupied cell's open interior.
  expectCollision(map, gridPoint(map, 2.5, 3.5), gridPoint(map, 4.5, 3.5));
  expectCollision(map, gridPoint(map, 3.0, 3.0), gridPoint(map, 4.0, 4.0));

  // Entering from the boundary, and connecting two opposite boundaries, are
  // collisions even though their endpoints themselves are legal contacts.
  expectCollision(map, gridPoint(map, 3.0, 3.5), gridPoint(map, 3.5, 3.5));
  expectCollision(map, gridPoint(map, 3.0, 3.5), gridPoint(map, 4.0, 3.5));

  // A genuine 1e-9-cell intrusion must not be hidden by numeric tolerances.
  expectCollision(map, gridPoint(map, 3.0 + 1.0e-9, 3.0), gridPoint(map, 3.0 + 1.0e-9, 4.0));

  // The occupied union is closed, but its external boundary is deliberately
  // legal for this planner's visibility segments.
  expectFree(map, gridPoint(map, 3.0, 3.0), gridPoint(map, 4.0, 3.0));
  expectFree(map, gridPoint(map, 3.0, 3.0), gridPoint(map, 3.0, 4.0));
  expectFree(map, gridPoint(map, 2.5, 3.5), gridPoint(map, 3.0, 3.5));
  expectFree(map, gridPoint(map, 3.0, 3.5), gridPoint(map, 2.5, 3.5));

  // This line touches only the lower-left occupied-cell vertex.
  expectFree(map, gridPoint(map, 2.0, 4.0), gridPoint(map, 4.0, 2.0));

  // Degenerate segments exercise point classification explicitly.
  expectCollision(map, gridPoint(map, 3.5, 3.5), gridPoint(map, 3.5, 3.5));
  expectFree(map, gridPoint(map, 3.0, 3.0), gridPoint(map, 3.0, 3.0));
}

TEST(PathCollisionOracle, CanonicalizesOnlyRepresentableBoundaryRoundTrips) {
  auto map = makeGrid(8, 8, 0.1F, 0.1, -0.2);
  setCell(map, 3, 3, 100);

  const double exact_boundary_y =
    map.info.origin.position.y + 3.0 * static_cast<double>(map.info.resolution);
  const double below = std::nextafter(exact_boundary_y, -std::numeric_limits<double>::infinity());
  const double above = std::nextafter(exact_boundary_y, std::numeric_limits<double>::infinity());

  auto start = gridPoint(map, 3.25, 3.0);
  auto end = gridPoint(map, 3.75, 3.0);
  start.y = below;
  end.y = below;
  expectFree(map, start, end);
  start.y = above;
  end.y = above;
  expectFree(map, start, end);

  // Moving well beyond the two-ULP round-trip window is a real intrusion.
  start.y = exact_boundary_y + 1.0e-9;
  end.y = exact_boundary_y + 1.0e-9;
  expectCollision(map, start, end);
}

TEST(PathCollisionOracle, RejectsOccupiedInternalHorizontalAndVerticalSeams) {
  auto vertical = makeGrid();
  setCell(vertical, 3, 3, 100);
  setCell(vertical, 4, 3, 100);
  expectCollision(vertical, gridPoint(vertical, 4.0, 3.2), gridPoint(vertical, 4.0, 3.8));
  // The seam's endpoint is incident to two free cells, so it is an external
  // boundary vertex rather than an interior point of the occupied union.
  expectFree(vertical, gridPoint(vertical, 4.0, 3.0), gridPoint(vertical, 4.0, 3.0));

  auto horizontal = makeGrid();
  setCell(horizontal, 3, 3, 100);
  setCell(horizontal, 3, 4, 100);
  expectCollision(horizontal, gridPoint(horizontal, 3.2, 4.0), gridPoint(horizontal, 3.8, 4.0));
}

TEST(PathCollisionOracle, FourOccupiedCellsMakeTheirSharedVertexInterior) {
  auto map = makeGrid();
  setCell(map, 3, 3, 100);
  setCell(map, 4, 3, 100);
  setCell(map, 3, 4, 100);
  setCell(map, 4, 4, 100);
  expectCollision(map, gridPoint(map, 4.0, 4.0), gridPoint(map, 4.0, 4.0));

  // With one free incident cell the same lattice vertex is on the external
  // boundary and therefore remains legal.
  setCell(map, 4, 4, 0);
  expectFree(map, gridPoint(map, 4.0, 4.0), gridPoint(map, 4.0, 4.0));
}

TEST(PathCollisionOracle, MatchesThresholdUnknownAndForcedOuterRingSemantics) {
  auto map = makeGrid();
  const auto start = gridPoint(map, 2.5, 3.5);
  const auto end = gridPoint(map, 4.5, 3.5);
  OracleOptions options;
  options.occupied_threshold = 50;

  setCell(map, 3, 3, 49);
  expectFree(map, start, end, options);
  setCell(map, 3, 3, 50);
  expectCollision(map, start, end, options);

  setCell(map, 3, 3, -1);
  options.allow_unknown = false;
  expectCollision(map, start, end, options);
  options.allow_unknown = true;
  expectFree(map, start, end, options);

  // Input values cannot make the planner's forced outer ring traversable.
  expectCollision(map, gridPoint(map, 0.5, 2.5), gridPoint(map, 0.5, 3.5), options);
  // Nor can a path leave the map: OOB has no free-cell closure.
  expectCollision(map, gridPoint(map, 1.5, 2.5), gridPoint(map, -0.5, 2.5), options);
}

TEST(PathCollisionOracle, UsesContinuousGridCoordinatesForMetricMaps) {
  auto map = makeGrid(12, 10, 0.25F, -4.0, 7.0);
  setCell(map, 5, 4, 100);

  expectCollision(map, gridPoint(map, 4.25, 4.5), gridPoint(map, 6.75, 4.5));
  expectFree(map, gridPoint(map, 5.125, 4.0), gridPoint(map, 5.875, 4.0));
}

TEST(PathCollisionOracle, ValidatesEveryPathSegmentAndReturnsDiagnostics) {
  auto map = makeGrid();
  setCell(map, 3, 3, 100);

  nav_msgs::msg::Path path;
  path.header.frame_id = map.header.frame_id;
  for (const auto& point : std::vector<geometry_msgs::msg::Point>{
         gridPoint(map, 2.0, 2.0), gridPoint(map, 2.5, 3.5), gridPoint(map, 4.5, 3.5)}) {
    geometry_msgs::msg::PoseStamped pose;
    pose.header.frame_id = map.header.frame_id;
    pose.pose.position = point;
    path.poses.push_back(pose);
  }

  const auto collision = validatePath(map, path);
  EXPECT_EQ(collision.status, ValidationStatus::kCollision);
  EXPECT_EQ(collision.segment_index, 1U);
  EXPECT_TRUE(std::isfinite(collision.grid_x));
  EXPECT_TRUE(std::isfinite(collision.grid_y));
  EXPECT_NE(collision.diagnostic.find("grid coordinate"), std::string::npos);

  path.poses.back().pose.position = gridPoint(map, 4.0, 3.0);
  path.poses[1].pose.position = gridPoint(map, 3.0, 3.0);
  const auto free = validatePath(map, path);
  EXPECT_EQ(free.status, ValidationStatus::kCollisionFree) << free.diagnostic;
}

TEST(PathCollisionOracle, RejectsMismatchedPathAndPoseFramesAsInconclusive) {
  auto map = makeGrid();
  nav_msgs::msg::Path path;
  path.header.frame_id = map.header.frame_id;
  for (double x : {2.0, 3.0}) {
    geometry_msgs::msg::PoseStamped pose;
    pose.header.frame_id = map.header.frame_id;
    pose.pose.position = gridPoint(map, x, 2.0);
    path.poses.push_back(pose);
  }

  path.header.frame_id = "odom";
  auto result = validatePath(map, path);
  EXPECT_EQ(result.status, ValidationStatus::kInconclusive);
  EXPECT_NE(result.diagnostic.find("Path frame_id"), std::string::npos);

  path.header.frame_id = map.header.frame_id;
  path.poses[1].header.frame_id = "odom";
  result = validatePath(map, path);
  EXPECT_EQ(result.status, ValidationStatus::kInconclusive);
  EXPECT_NE(result.diagnostic.find("Path pose 1"), std::string::npos);

  path.poses[1].header.frame_id = map.header.frame_id;
  map.header.frame_id.clear();
  result = validatePath(map, path);
  EXPECT_EQ(result.status, ValidationStatus::kInconclusive);
  EXPECT_NE(result.diagnostic.find("map frame_id"), std::string::npos);
}

TEST(PathCollisionOracle, MatchesProductionIdentityQuaternionTolerance) {
  auto map = makeGrid();
  map.info.origin.orientation.w = 1.0 + 1.0e-7;
  const auto point = gridPoint(map, 2.0, 2.0);
  const auto result = validateSegment(map, point, point);
  EXPECT_EQ(result.status, ValidationStatus::kInconclusive);
  EXPECT_NE(result.diagnostic.find("identity map rotation"), std::string::npos);
}

TEST(PathCollisionOracle, ReportsUnrepresentableMetricGridAsInconclusive) {
  auto map = makeGrid(8, 8, 0.25F, 1.0e16, 1.0e16);
  const auto point = gridPoint(map, 2.0, 2.0);
  auto result = validateSegment(map, point, point);
  EXPECT_EQ(result.status, ValidationStatus::kInconclusive);
  EXPECT_FALSE(result.conclusive());
  EXPECT_NE(result.diagnostic.find("not distinguishable"), std::string::npos);

  // Distinct boundary doubles are still insufficient when their two-ULP
  // canonicalization windows overlap.  At 2^53, a four-meter grid has a
  // two-meter ULP: origin + 14 represents grid x=3.5 but is only one ULP from
  // the x=4 boundary.  The oracle must decline this ambiguous map instead of
  // snapping a genuine occupied-cell interior point onto a legal boundary.
  map = makeGrid(8, 8, 4.0F, 9007199254740992.0, 0.0);
  setCell(map, 3, 3, 100);
  auto start = gridPoint(map, 3.5, 3.25);
  auto end = gridPoint(map, 3.5, 3.75);
  result = validateSegment(map, start, end);
  EXPECT_EQ(result.status, ValidationStatus::kInconclusive);
  EXPECT_FALSE(result.conclusive());
  EXPECT_NE(result.diagnostic.find("two-ULP snap windows"), std::string::npos);
}

TEST(PathSelfIntersectionOracle, AcceptsSimpleInterpolatedPolylineAndAdjacentSharedEndpoints) {
  const auto simple =
    makePath({{0.0, 0.0}, {0.5, 0.0}, {1.0, 0.0}, {1.5, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {3.0, 1.0}});
  const auto result = validateNoSelfIntersection(simple);
  EXPECT_EQ(result.status, SelfIntersectionStatus::kIntersectionFree) << result.diagnostic;
  EXPECT_EQ(result.kind, SelfIntersectionKind::kNone);
  EXPECT_TRUE(result.intersectionFree());
  EXPECT_TRUE(result.conclusive());
  EXPECT_FALSE(result.diagnostic.empty());

  // A single degenerate segment has no non-consecutive segment to intersect.
  const auto degenerate = validateNoSelfIntersection(makePath({{1.0, 1.0}, {1.0, 1.0}}));
  EXPECT_EQ(degenerate.status, SelfIntersectionStatus::kIntersectionFree) << degenerate.diagnostic;
}

TEST(PathSelfIntersectionOracle, ClassifiesProperCrossingExactly) {
  const auto result =
    validateNoSelfIntersection(makePath({{0.0, 0.0}, {4.0, 4.0}, {0.0, 4.0}, {4.0, 0.0}}));
  EXPECT_EQ(result.status, SelfIntersectionStatus::kSelfIntersection) << result.diagnostic;
  EXPECT_EQ(result.kind, SelfIntersectionKind::kProperCrossing);
  EXPECT_EQ(result.first_segment_index, 0U);
  EXPECT_EQ(result.second_segment_index, 2U);
  EXPECT_TRUE(result.selfIntersectionDetected());
  EXPECT_FALSE(result.diagnostic.empty());
}

TEST(PathSelfIntersectionOracle, ClassifiesNonAdjacentInteriorTouch) {
  const auto result =
    validateNoSelfIntersection(makePath({{0.0, 0.0}, {4.0, 0.0}, {4.0, 4.0}, {2.0, 0.0}}));
  EXPECT_EQ(result.status, SelfIntersectionStatus::kSelfIntersection) << result.diagnostic;
  EXPECT_EQ(result.kind, SelfIntersectionKind::kNonAdjacentTouch);
  EXPECT_EQ(result.first_segment_index, 0U);
  EXPECT_EQ(result.second_segment_index, 2U);
}

TEST(PathSelfIntersectionOracle, ClassifiesRepeatedNonConsecutiveVertex) {
  const auto result = validateNoSelfIntersection(
    makePath({{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {2.0, 0.0}, {4.0, 0.0}}));
  EXPECT_EQ(result.status, SelfIntersectionStatus::kSelfIntersection) << result.diagnostic;
  EXPECT_EQ(result.kind, SelfIntersectionKind::kRepeatedVertex);
  EXPECT_EQ(result.first_segment_index, 0U);
  EXPECT_EQ(result.second_segment_index, 2U);
}

TEST(PathSelfIntersectionOracle, ClassifiesPositiveLengthCollinearOverlap) {
  const auto result =
    validateNoSelfIntersection(makePath({{0.0, 0.0}, {4.0, 0.0}, {3.0, 0.0}, {1.0, 0.0}}));
  EXPECT_EQ(result.status, SelfIntersectionStatus::kSelfIntersection) << result.diagnostic;
  EXPECT_EQ(result.kind, SelfIntersectionKind::kCollinearOverlap);
  EXPECT_EQ(result.first_segment_index, 0U);
  EXPECT_EQ(result.second_segment_index, 2U);
}

TEST(PathSelfIntersectionOracle, TreatsClosedPathContactAsRepeatedVertex) {
  const auto result = validateNoSelfIntersection(
    makePath({{0.0, 0.0}, {2.0, 0.0}, {2.0, 2.0}, {0.0, 2.0}, {0.0, 0.0}}));
  EXPECT_EQ(result.status, SelfIntersectionStatus::kSelfIntersection) << result.diagnostic;
  EXPECT_EQ(result.kind, SelfIntersectionKind::kRepeatedVertex);
  EXPECT_EQ(result.first_segment_index, 0U);
  EXPECT_EQ(result.second_segment_index, 3U);
}

TEST(PathSelfIntersectionOracle, RejectsAmbiguousInputInsteadOfGuessing) {
  auto path = makePath({{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}});
  path.poses[1].pose.position.x = std::numeric_limits<double>::quiet_NaN();
  auto result = validateNoSelfIntersection(path);
  EXPECT_EQ(result.status, SelfIntersectionStatus::kInconclusive);
  EXPECT_FALSE(result.conclusive());
  EXPECT_NE(result.diagnostic.find("non-finite"), std::string::npos);

  path = makePath({{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}});
  path.poses[1].header.frame_id = "odom";
  result = validateNoSelfIntersection(path);
  EXPECT_EQ(result.status, SelfIntersectionStatus::kInconclusive);
  EXPECT_NE(result.diagnostic.find("frame_id"), std::string::npos);

  path = makePath({{0.0, 0.0}, {1.0, 1.0}, {2.0, 0.0}});
  result = validateNoSelfIntersection(path, 1);
  EXPECT_EQ(result.status, SelfIntersectionStatus::kInconclusive);
  EXPECT_NE(result.diagnostic.find("max_segments"), std::string::npos);
}

}  // namespace
