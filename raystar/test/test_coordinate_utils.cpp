#include <gtest/gtest.h>
#include <raystar/coordinate_utils.h>

#include "metric_bound_search.h"

#include <cfenv>
#include <cmath>
#include <limits>

using namespace raystar;

namespace {

void expectMinimalGridSearchSuperset(double metric_bound, float resolution) {
  double grid_bound = std::numeric_limits<double>::quiet_NaN();
  ASSERT_TRUE(gridCostSearchUpperBoundForMetricBound(metric_bound, resolution, grid_bound));
  if (detail::compareExactMetricProduct(
        std::numeric_limits<double>::max(), resolution, metric_bound) >= 0) {
    EXPECT_GE(detail::compareExactMetricProduct(grid_bound, resolution, metric_bound), 0);
    if (grid_bound > 0.0) {
      const double predecessor = std::nextafter(grid_bound, 0.0);
      EXPECT_LT(detail::compareExactMetricProduct(predecessor, resolution, metric_bound), 0);
    }
  } else {
    EXPECT_DOUBLE_EQ(grid_bound, std::numeric_limits<double>::max());
  }
}

}  // namespace

TEST(CoordinateUtils, PublishedCostBoundCalibrationHandlesBothDoubleRoundingDirections) {
  const float resolution = 0.1f;

  const double reusable_grid_cost = std::hypot(6.0, 5.0);
  const double reusable_metric_bound = gridCostToMetric(reusable_grid_cost, resolution);
  double reusable_grid_bound = 0.0;
  ASSERT_TRUE(
    gridCostSearchUpperBoundForMetricBound(reusable_metric_bound, resolution, reusable_grid_bound));
  EXPECT_GE(reusable_grid_bound, reusable_grid_cost);
  expectMinimalGridSearchSuperset(reusable_metric_bound, resolution);

  const double excluded_grid_cost = std::hypot(4.0, 1.0);
  const double excluded_metric_bound =
    std::nextafter(gridCostToMetric(excluded_grid_cost, resolution), 0.0);
  double excluded_grid_bound = 0.0;
  ASSERT_TRUE(
    gridCostSearchUpperBoundForMetricBound(excluded_metric_bound, resolution, excluded_grid_bound));
  EXPECT_LE(excluded_grid_bound, excluded_grid_cost);
  expectMinimalGridSearchSuperset(excluded_metric_bound, resolution);

  // The tight metric ceiling of this direct path is one ULP below the rounded
  // fma of its tight grid ceiling. A publication-limited conversion rejects
  // the path, while the exact quotient ceiling correctly keeps it searchable.
  constexpr float fine_resolution = 0x1.99999ap-5f;
  constexpr double metric_path_bound = 0x1.cd82b4b9764c5p-4;
  constexpr double grid_path_ceiling = 0x1.2071b0abcd839p+1;
  double search_bound = 0.0;
  ASSERT_TRUE(
    gridCostSearchUpperBoundForMetricBound(metric_path_bound, fine_resolution, search_bound));
  EXPECT_DOUBLE_EQ(search_bound, grid_path_ceiling);
  EXPECT_GT(gridCostToMetric(grid_path_ceiling, fine_resolution), metric_path_bound);
  expectMinimalGridSearchSuperset(metric_path_bound, fine_resolution);
}

TEST(CoordinateUtils, PublishedCostBoundCalibrationHandlesSubnormalAndOverflowEdges) {
  const double minimum_metric = std::numeric_limits<double>::denorm_min();
  const float minimum_resolution = std::numeric_limits<float>::denorm_min();
  expectMinimalGridSearchSuperset(minimum_metric, minimum_resolution);

  double zero_grid_bound = 1.0;
  ASSERT_TRUE(gridCostSearchUpperBoundForMetricBound(
    minimum_metric, std::numeric_limits<float>::max(), zero_grid_bound));
  EXPECT_DOUBLE_EQ(zero_grid_bound, std::numeric_limits<double>::denorm_min());
  expectMinimalGridSearchSuperset(minimum_metric, std::numeric_limits<float>::max());

  double saturated_grid_bound = 0.0;
  ASSERT_TRUE(gridCostSearchUpperBoundForMetricBound(
    std::numeric_limits<double>::max(), minimum_resolution, saturated_grid_bound));
  EXPECT_DOUBLE_EQ(saturated_grid_bound, std::numeric_limits<double>::max());
  expectMinimalGridSearchSuperset(std::numeric_limits<double>::max(), minimum_resolution);

  const float two = 2.0f;
  double overflow_limited_grid_bound = 0.0;
  ASSERT_TRUE(gridCostSearchUpperBoundForMetricBound(
    std::numeric_limits<double>::max(), two, overflow_limited_grid_bound));
  EXPECT_EQ(std::ilogb(overflow_limited_grid_bound), 1022);
  expectMinimalGridSearchSuperset(std::numeric_limits<double>::max(), two);
}

TEST(CoordinateUtils, PublishedCostBoundCalibrationRejectsInvalidInputs) {
  double grid_bound = 123.0;
  EXPECT_FALSE(gridCostSearchUpperBoundForMetricBound(
    std::numeric_limits<double>::quiet_NaN(), 1.0f, grid_bound));
  EXPECT_TRUE(std::isnan(grid_bound));
  EXPECT_FALSE(gridCostSearchUpperBoundForMetricBound(
    std::numeric_limits<double>::infinity(), 1.0f, grid_bound));
  EXPECT_FALSE(gridCostSearchUpperBoundForMetricBound(-1.0, 1.0f, grid_bound));
  EXPECT_FALSE(gridCostSearchUpperBoundForMetricBound(1.0, 0.0f, grid_bound));
  EXPECT_FALSE(gridCostSearchUpperBoundForMetricBound(
    1.0, std::numeric_limits<float>::infinity(), grid_bound));
}

TEST(CoordinateUtils, MetricSearchPaddingCoversMaxNodesFrontierChildInEveryRoundingMode) {
  const int original_rounding = std::fegetround();
  ASSERT_NE(original_rounding, -1);
  struct RoundingGuard {
    int mode;
    ~RoundingGuard() {
      (void)std::fesetround(mode);
    }
  } rounding_guard{original_rounding};
  (void)rounding_guard;

  GridMap map;
  map.width = 2;
  map.height = 2;
  map.resolution = 1.0f;
  // The world interval crosses zero, so endpoint FMA errors can have opposite
  // signs and the full two-endpoint contraction proof is required.
  volatile double runtime_origin = -1.0;
  map.origin_x = runtime_origin;
  map.origin_y = runtime_origin;

  constexpr double kSquareRootTwoUpper = 0x1.6a09e667f3bcdp+0;
  const long double operation_spacing = std::ldexp(1.0L, -52);
  // max_nodes=1 still permits a queued child optimistic certificate with two
  // segments: start->child and child->goal.
  const long double required_two_segment_contraction =
    4.0L * static_cast<long double>(kSquareRootTwoUpper) * operation_spacing;

  for (const int rounding_mode : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
    SCOPED_TRACE(rounding_mode);
    ASSERT_EQ(std::fesetround(rounding_mode), 0);
    volatile double runtime_bound = 1.0;
    double padded_bound = 0.0;
    ASSERT_TRUE(paddedMetricBoundForGridSearch(map, runtime_bound, 1u, padded_bound));
    EXPECT_TRUE(std::isfinite(padded_bound));
    EXPECT_GE(static_cast<long double>(padded_bound) - 1.0L, required_two_segment_contraction);
  }
}

TEST(CoordinateUtils, MetricSearchPaddingSaturatesAtMaximumInEveryRoundingMode) {
  const int original_rounding = std::fegetround();
  ASSERT_NE(original_rounding, -1);
  struct RoundingGuard {
    int mode;
    ~RoundingGuard() {
      (void)std::fesetround(mode);
    }
  } rounding_guard{original_rounding};
  (void)rounding_guard;

  GridMap map;
  map.width = 2;
  map.height = 2;
  map.resolution = 1.0f;
  map.origin_x = -1.0;
  map.origin_y = -1.0;
  volatile double runtime_maximum = std::numeric_limits<double>::max();
  for (const int rounding_mode : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
    SCOPED_TRACE(rounding_mode);
    ASSERT_EQ(std::fesetround(rounding_mode), 0);
    double padded_bound = 0.0;
    ASSERT_TRUE(paddedMetricBoundForGridSearch(map, runtime_maximum, 1u, padded_bound));
    EXPECT_DOUBLE_EQ(padded_bound, std::numeric_limits<double>::max());
  }
}

TEST(CoordinateUtils, WorldToMapBasic) {
  GridMap map;
  map.width = 100;
  map.height = 100;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.resize(100 * 100, 0);

  unsigned int mx, my;
  ASSERT_TRUE(worldToMap(map, 50.0, 50.0, mx, my));
  EXPECT_EQ(mx, 50u);
  EXPECT_EQ(my, 50u);

  ASSERT_TRUE(worldToMap(map, 0.0, 0.0, mx, my));
  EXPECT_EQ(mx, 0u);
  EXPECT_EQ(my, 0u);

  ASSERT_TRUE(worldToMap(map, 99.9, 99.9, mx, my));
  EXPECT_EQ(mx, 99u);
  EXPECT_EQ(my, 99u);
}

TEST(CoordinateUtils, WorldToMapOutOfRange) {
  GridMap map;
  map.width = 10;
  map.height = 10;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.resize(100, 0);

  unsigned int mx, my;
  EXPECT_FALSE(worldToMap(map, -1.0, 0.0, mx, my));
  EXPECT_FALSE(worldToMap(map, 0.0, -1.0, mx, my));
  EXPECT_FALSE(worldToMap(map, 10.0, 0.0, mx, my));
  EXPECT_FALSE(worldToMap(map, 0.0, 10.0, mx, my));
}

TEST(CoordinateUtils, WorldToMapRejectsInvalidMetricMetadataBeforeCasting) {
  GridMap valid_map;
  valid_map.width = 10;
  valid_map.height = 10;
  valid_map.resolution = 1.0f;
  valid_map.origin_x = 0.0;
  valid_map.origin_y = 0.0;
  valid_map.data.resize(100, 0);

  const auto expect_rejected = [](const GridMap& map, double wx, double wy) {
    unsigned int mx = 123u;
    unsigned int my = 456u;
    EXPECT_FALSE(worldToMap(map, wx, wy, mx, my));
    EXPECT_EQ(mx, 123u) << "A rejected conversion must not modify mx";
    EXPECT_EQ(my, 456u) << "A rejected conversion must not modify my";
  };

  auto invalid_map = valid_map;
  invalid_map.resolution = 0.0f;
  expect_rejected(invalid_map, 1.0, 1.0);

  invalid_map = valid_map;
  invalid_map.resolution = -1.0f;
  expect_rejected(invalid_map, -1.0, -1.0);

  invalid_map = valid_map;
  invalid_map.resolution = std::numeric_limits<float>::quiet_NaN();
  expect_rejected(invalid_map, 1.0, 1.0);

  invalid_map = valid_map;
  invalid_map.resolution = std::numeric_limits<float>::infinity();
  expect_rejected(invalid_map, 1.0, 1.0);

  invalid_map = valid_map;
  invalid_map.origin_x = std::numeric_limits<double>::quiet_NaN();
  expect_rejected(invalid_map, 1.0, 1.0);

  invalid_map = valid_map;
  invalid_map.origin_y = std::numeric_limits<double>::infinity();
  expect_rejected(invalid_map, 1.0, 1.0);

  invalid_map = valid_map;
  invalid_map.width = 0;
  expect_rejected(invalid_map, 0.0, 0.0);

  invalid_map = valid_map;
  invalid_map.height = 0;
  expect_rejected(invalid_map, 0.0, 0.0);
}

TEST(CoordinateUtils, WorldToMapRejectsNonFiniteAndUnrepresentableCoordinates) {
  GridMap map;
  map.width = 10;
  map.height = 10;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.resize(100, 0);

  const auto expect_rejected = [&map](double wx, double wy) {
    unsigned int mx = 123u;
    unsigned int my = 456u;
    EXPECT_FALSE(worldToMap(map, wx, wy, mx, my));
    EXPECT_EQ(mx, 123u) << "A rejected conversion must not modify mx";
    EXPECT_EQ(my, 456u) << "A rejected conversion must not modify my";
  };

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  const double maximum = std::numeric_limits<double>::max();
  expect_rejected(nan, 1.0);
  expect_rejected(1.0, nan);
  expect_rejected(infinity, 1.0);
  expect_rejected(1.0, infinity);
  expect_rejected(-infinity, 1.0);
  expect_rejected(1.0, -infinity);
  expect_rejected(maximum, 1.0);
  expect_rejected(1.0, maximum);

  auto overflow_map = map;
  overflow_map.origin_x = -maximum;
  unsigned int mx = 123u;
  unsigned int my = 456u;
  EXPECT_FALSE(worldToMap(overflow_map, maximum, 1.0, mx, my));
  EXPECT_EQ(mx, 123u);
  EXPECT_EQ(my, 456u);
}

TEST(CoordinateUtils, WorldToMapKeepsValidFractionalResolutionSemantics) {
  GridMap map;
  map.width = 4;
  map.height = 5;
  map.resolution = 0.5f;
  map.origin_x = -1.0;
  map.origin_y = 2.0;
  map.data.resize(20, 0);

  unsigned int mx = 0;
  unsigned int my = 0;
  ASSERT_TRUE(worldToMap(map, -0.01, 2.99, mx, my));
  EXPECT_EQ(mx, 1u);
  EXPECT_EQ(my, 1u);

  EXPECT_FALSE(worldToMap(map, 1.0, 2.0, mx, my));
  EXPECT_FALSE(worldToMap(map, -1.0, 4.5, mx, my));
}

TEST(CoordinateUtils, MapToWorldAndBack) {
  GridMap map;
  map.width = 100;
  map.height = 100;
  map.resolution = 1.0f;
  map.origin_x = -50.0;
  map.origin_y = -50.0;
  map.data.resize(100 * 100, 0);

  for (unsigned int mx = 0; mx < 100; mx += 10) {
    for (unsigned int my = 0; my < 100; my += 10) {
      double wx, wy;
      ASSERT_TRUE(mapToWorld(map, mx, my, wx, wy));
      unsigned int mx2, my2;
      ASSERT_TRUE(worldToMap(map, wx, wy, mx2, my2)) << "mx=" << mx << " my=" << my;
      EXPECT_EQ(mx, mx2) << "wx=" << wx;
      EXPECT_EQ(my, my2) << "wy=" << wy;
    }
  }
}

TEST(CoordinateUtils, WorldToContinuousMapPreservesFractionAndFloorCell) {
  GridMap map;
  map.width = 20;
  map.height = 20;
  map.resolution = 0.25f;
  map.origin_x = -4.0;
  map.origin_y = 7.0;
  map.data.assign(400, 0);

  ContinuousGridPoint point;
  ASSERT_TRUE(worldToContinuousMap(map, -3.4375, 7.875, point));
  EXPECT_DOUBLE_EQ(point.x, 2.25);
  EXPECT_DOUBLE_EQ(point.y, 3.5);
  EXPECT_EQ(point.cell_x, 2u);
  EXPECT_EQ(point.cell_y, 3u);

  double wx = 0.0;
  double wy = 0.0;
  ASSERT_TRUE(continuousMapToWorld(map, point.x, point.y, wx, wy));
  EXPECT_DOUBLE_EQ(wx, -3.4375);
  EXPECT_DOUBLE_EQ(wy, 7.875);

  const ContinuousGridPoint original = point;
  EXPECT_FALSE(worldToContinuousMap(map, -4.0, 12.0, point));
  EXPECT_DOUBLE_EQ(point.x, original.x);
  EXPECT_DOUBLE_EQ(point.y, original.y);
  EXPECT_EQ(point.cell_x, original.cell_x);
  EXPECT_EQ(point.cell_y, original.cell_y);
}

TEST(CoordinateUtils, SnapsRepresentableWorldGridBoundariesBeforeFlooring) {
  GridMap map;
  map.width = 10;
  map.height = 10;
  // 0.1 is intentionally not exactly representable as a float.  Construct
  // the world point through the same forward transform used by the node.
  map.resolution = 0.1f;
  map.origin_x = 0.1;
  map.origin_y = -0.2;
  map.data.assign(100, 0);

  const double boundary_x = map.origin_x + 2.0 * static_cast<double>(map.resolution);
  const double boundary_y = map.origin_y + 3.0 * static_cast<double>(map.resolution);

  ContinuousGridPoint point;
  ASSERT_TRUE(worldToContinuousMap(map, boundary_x, boundary_y, point));
  EXPECT_DOUBLE_EQ(point.x, 2.0);
  EXPECT_DOUBLE_EQ(point.y, 3.0);
  EXPECT_EQ(point.cell_x, 2u);
  EXPECT_EQ(point.cell_y, 3u);

  unsigned int cell_x = 0;
  unsigned int cell_y = 0;
  ASSERT_TRUE(worldToMap(map, boundary_x, boundary_y, cell_x, cell_y));
  EXPECT_EQ(cell_x, 2u);
  EXPECT_EQ(cell_y, 3u);
}

TEST(CoordinateUtils, RejectsRepresentableWorldUpperBoundaryAfterSnapping) {
  GridMap map;
  map.width = 10;
  map.height = 10;
  map.resolution = 0.1f;
  map.origin_x = 0.1;
  map.origin_y = -0.2;
  map.data.assign(100, 0);

  const double map_right =
    map.origin_x + static_cast<double>(map.width) * static_cast<double>(map.resolution);
  const double map_top =
    map.origin_y + static_cast<double>(map.height) * static_cast<double>(map.resolution);
  ContinuousGridPoint point{7.0, 8.0, 7u, 8u};
  EXPECT_FALSE(worldToContinuousMap(map, map_right, 0.0, point));
  EXPECT_DOUBLE_EQ(point.x, 7.0);
  EXPECT_DOUBLE_EQ(point.y, 8.0);
  EXPECT_EQ(point.cell_x, 7u);
  EXPECT_EQ(point.cell_y, 8u);
  EXPECT_FALSE(worldToContinuousMap(map, 0.0, map_top, point));
}

TEST(CoordinateUtils, MapToWorldRejectsInvalidTransformAndReturnsNaN) {
  GridMap valid_map;
  valid_map.width = 10;
  valid_map.height = 10;
  valid_map.resolution = 1.0f;
  valid_map.origin_x = 0.0;
  valid_map.origin_y = 0.0;
  valid_map.data.resize(100, 0);

  const auto expect_rejected = [](const GridMap& map) {
    double wx = 123.0;
    double wy = 456.0;
    EXPECT_FALSE(mapToWorld(map, 1u, 2u, wx, wy));
    EXPECT_TRUE(std::isnan(wx));
    EXPECT_TRUE(std::isnan(wy));
  };

  auto invalid_map = valid_map;
  invalid_map.resolution = 0.0f;
  expect_rejected(invalid_map);

  invalid_map = valid_map;
  invalid_map.resolution = -1.0f;
  expect_rejected(invalid_map);

  invalid_map = valid_map;
  invalid_map.resolution = std::numeric_limits<float>::quiet_NaN();
  expect_rejected(invalid_map);

  invalid_map = valid_map;
  invalid_map.resolution = std::numeric_limits<float>::infinity();
  expect_rejected(invalid_map);

  invalid_map = valid_map;
  invalid_map.origin_x = std::numeric_limits<double>::quiet_NaN();
  expect_rejected(invalid_map);

  invalid_map = valid_map;
  invalid_map.origin_y = std::numeric_limits<double>::infinity();
  expect_rejected(invalid_map);
}

TEST(CoordinateUtils, GridMapAccess) {
  GridMap map;
  map.width = 10;
  map.height = 10;
  map.data.resize(100, 0);

  map(5, 5) = 1;
  EXPECT_EQ(map(5, 5), 1);
  EXPECT_EQ(map(0, 0), 0);
  EXPECT_EQ(map(9u, 9u), 0);
}

TEST(CoordinateUtils, GridMapOutOfBoundsReturnsOccupied) {
  GridMap map;
  map.width = 10;
  map.height = 10;
  map.data.resize(100, 0);

  EXPECT_EQ(map.at(10u, 0u), 1u);
  EXPECT_EQ(map.at(0u, 10u), 1u);
}

TEST(CoordinateUtils, NormalizeAngle) {
  EXPECT_NEAR(normalize_angle(0.0), 0.0, 1e-10);
  EXPECT_NEAR(normalize_angle(kPi), kPi, 1e-10);
  EXPECT_NEAR(normalize_angle(-kPi), kPi, 1e-10);
  EXPECT_NEAR(normalize_angle(kTwoPi), 0.0, 1e-10);
  EXPECT_NEAR(normalize_angle(-kTwoPi), 0.0, 1e-10);
  EXPECT_NEAR(normalize_angle(kPi / 2.0), kPi / 2.0, 1e-10);
  EXPECT_NEAR(normalize_angle(-kPi / 2.0), -kPi / 2.0, 1e-10);
}

TEST(CoordinateUtils, NormalizeAnglePositive) {
  EXPECT_NEAR(normalize_angle_positive(0.0), 0.0, 1e-10);
  EXPECT_NEAR(normalize_angle_positive(kPi), kPi, 1e-10);
  EXPECT_NEAR(normalize_angle_positive(-kPi), kPi, 1e-10);
  EXPECT_NEAR(normalize_angle_positive(kTwoPi), 0.0, 1e-10);
  EXPECT_NEAR(normalize_angle_positive(kPi / 2.0), kPi / 2.0, 1e-10);
}

TEST(CoordinateUtils, AngleNormalizationHandlesNonFiniteAndExtremeFiniteInputs) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  for (const double non_finite : {nan, infinity, -infinity}) {
    EXPECT_TRUE(std::isnan(normalize_angle(non_finite)));
    EXPECT_TRUE(std::isnan(normalize_angle_positive(non_finite)));
  }

  for (const double extreme :
       {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()}) {
    const double signed_angle = normalize_angle(extreme);
    EXPECT_TRUE(std::isfinite(signed_angle));
    EXPECT_GT(signed_angle, -kPi);
    EXPECT_LE(signed_angle, kPi);

    const double positive_angle = normalize_angle_positive(extreme);
    EXPECT_TRUE(std::isfinite(positive_angle));
    EXPECT_GE(positive_angle, 0.0);
    EXPECT_LT(positive_angle, kTwoPi);
  }
}

TEST(CoordinateUtils, CounterClockwiseSweepHandlesAllSpanClasses) {
  // Span < PI.
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(0.0, 0.0, kPi / 2.0));
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(kPi / 4.0, 0.0, kPi / 2.0));
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(kPi / 2.0, 0.0, kPi / 2.0));
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(-5e-8, 0.0, kPi / 2.0));
  EXPECT_FALSE(isAngleWithinCounterClockwiseSweep(-2e-7, 0.0, kPi / 2.0));
  EXPECT_FALSE(isAngleWithinCounterClockwiseSweep(kPi, 0.0, kPi / 2.0));

  // Span == PI.
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(-kPi / 2.0, -kPi / 2.0, kPi / 2.0));
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(0.0, -kPi / 2.0, kPi / 2.0));
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(kPi / 2.0, -kPi / 2.0, kPi / 2.0));
  EXPECT_FALSE(isAngleWithinCounterClockwiseSweep(-3.0 * kPi / 4.0, -kPi / 2.0, kPi / 2.0));

  // Span > PI: this is the H-01 regression that signed normalization clipped.
  const double wide_start = 0.25;
  const double wide_end = wide_start + 1.5 * kPi;
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(wide_start, wide_start, wide_end));
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(wide_start + kPi + 0.1, wide_start, wide_end));
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(wide_end, wide_start, wide_end));
  EXPECT_FALSE(isAngleWithinCounterClockwiseSweep(wide_start - 0.1, wide_start, wide_end));

  // Sweep crosses the -PI/PI branch cut.
  const double wrap_start = 3.0 * kPi / 4.0;
  const double wrap_end = -3.0 * kPi / 4.0;
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(wrap_start, wrap_start, wrap_end));
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(wrap_end, wrap_start, wrap_end));
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(kPi, wrap_start, wrap_end));
  EXPECT_TRUE(isAngleWithinCounterClockwiseSweep(-7.0 * kPi / 8.0, wrap_start, wrap_end));
  EXPECT_FALSE(isAngleWithinCounterClockwiseSweep(0.0, wrap_start, wrap_end));
}

TEST(SegmentsIntersect, ProperCrossing) {
  EXPECT_TRUE(segmentsIntersect({0, 0}, {2, 2}, {0, 2}, {2, 0}));
}

TEST(SegmentsIntersect, ParallelNoCrossing) {
  EXPECT_FALSE(segmentsIntersect({0, 0}, {1, 0}, {0, 1}, {1, 1}));
}

TEST(SegmentsIntersect, SharedEndpointTouches) {
  EXPECT_TRUE(segmentsIntersect({0, 0}, {1, 0}, {1, 0}, {2, 0}));
}

TEST(SegmentsIntersect, CollinearOverlap) {
  EXPECT_TRUE(segmentsIntersect({0, 0}, {2, 0}, {1, 0}, {3, 0}));
}

TEST(SegmentsIntersect, CollinearNoOverlap) {
  EXPECT_FALSE(segmentsIntersect({0, 0}, {1, 0}, {2, 0}, {3, 0}));
}

TEST(SegmentsIntersect, Backtracking) {
  EXPECT_TRUE(segmentsIntersect({0, 0}, {2, 0}, {2, 0}, {0, 0}));
}

TEST(SegmentsIntersect, EndpointOnInterior) {
  EXPECT_TRUE(segmentsIntersect({0, 0}, {2, 0}, {1, 0}, {1, 2}));
}

TEST(SegmentsIntersect, FullIntDomainDoesNotOverflowOrientationPredicates) {
  constexpr int lowest = std::numeric_limits<int>::lowest();
  constexpr int highest = std::numeric_limits<int>::max();

  EXPECT_EQ(orientation({lowest, lowest}, {highest, lowest}, {highest, highest}), 2);
  EXPECT_EQ(orientation({lowest, lowest}, {lowest, highest}, {highest, highest}), 1);
  EXPECT_EQ(orientation({lowest, lowest}, {0, 0}, {highest, highest}), 0);
  // Both determinant terms are non-zero and have opposite signs; their
  // mathematical difference is wider than a signed 64-bit integer.
  EXPECT_EQ(orientation({lowest, 0}, {0, highest}, {highest, lowest}), 1);
  EXPECT_EQ(orientation({lowest, 0}, {0, lowest}, {highest, highest}), 2);

  EXPECT_TRUE(
    segmentsIntersect({lowest, lowest}, {highest, highest}, {lowest, highest}, {highest, lowest}));
}

TEST(SegmentsIntersect, SharedVertexAtEnd) {
  EXPECT_TRUE(segmentsIntersect({0, 0}, {1, 1}, {2, 0}, {1, 1}));
}

TEST(NewPathSelfCrosses, ShortPathNoCrossing) {
  std::vector<std::pair<int, int>> path = {{0, 0}};
  EXPECT_FALSE(newPathSelfCrosses(path, {1, 1}));
}

TEST(NewPathSelfCrosses, NormalExtension) {
  std::vector<std::pair<int, int>> path = {{0, 0}, {2, 0}, {2, 2}};
  EXPECT_FALSE(newPathSelfCrosses(path, {0, 2}));
}

TEST(NewPathSelfCrosses, VertexRevisitStart) {
  std::vector<std::pair<int, int>> path = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
  EXPECT_TRUE(newPathSelfCrosses(path, {0, 0}));
}

TEST(NewPathSelfCrosses, VertexRevisitMiddle) {
  std::vector<std::pair<int, int>> path = {{0, 0}, {2, 0}, {2, 2}};
  EXPECT_TRUE(newPathSelfCrosses(path, {2, 0}));
}

TEST(NewPathSelfCrosses, GeometricCrossing) {
  std::vector<std::pair<int, int>> path = {{0, 0}, {2, 0}, {2, 2}, {0, 2}};
  EXPECT_TRUE(newPathSelfCrosses(path, {1, -1}));
}

TEST(NewPathSelfCrosses, ConsecutiveEndpointNotFlagged) {
  std::vector<std::pair<int, int>> path = {{0, 0}, {1, 0}};
  EXPECT_FALSE(newPathSelfCrosses(path, {2, 0}));
}

TEST(NewPathSelfCrosses, CooperativeOverloadPreservesLegacyResults) {
  const StopToken never_stop;

  bool self_crosses = true;
  EXPECT_EQ(newPathSelfCrosses({{0, 0}, {2, 0}, {2, 2}}, {0, 2}, self_crosses, never_stop),
            OperationStatus::success);
  EXPECT_FALSE(self_crosses);

  self_crosses = false;
  EXPECT_EQ(newPathSelfCrosses({{0, 0}, {2, 0}, {2, 2}}, {2, 0}, self_crosses, never_stop),
            OperationStatus::success);
  EXPECT_TRUE(self_crosses);
}

TEST(NewPathSelfCrosses, CooperativeStopDoesNotCommitPartialResult) {
  std::vector<std::pair<int, int>> path;
  for (int x = 0; x < 256; ++x) path.emplace_back(x, 0);

  size_t poll_count = 0;
  const StopToken stop_token([&poll_count]() { return ++poll_count == 8; });
  bool self_crosses = true;

  EXPECT_EQ(newPathSelfCrosses(path, {512, 1}, self_crosses, stop_token), OperationStatus::stopped);
  EXPECT_TRUE(self_crosses);
  EXPECT_EQ(poll_count, 8u);
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
