#include <gtest/gtest.h>

#include "small_grid_property.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <set>

namespace {

using raystar::test_property::describeSmallGridCase;
using raystar::test_property::fourNeighborReachable;
using raystar::test_property::GridCell;
using raystar::test_property::hasUnsupportedDiagonalContact;
using raystar::test_property::hasValidOccupancyValues;
using raystar::test_property::isEffectivelyBlocked;
using raystar::test_property::isForcedOuterRing;
using raystar::test_property::kOccupiedThreshold;
using raystar::test_property::kPropertyBaseSeed;
using raystar::test_property::kSeedStride;
using raystar::test_property::makeSmallGridCase;
using raystar::test_property::propertySeed;
using raystar::test_property::SmallGridCase;

constexpr std::size_t kCorpusSize = 24;

SmallGridCase makeUniformCase(std::uint32_t width,
                              std::uint32_t height,
                              GridCell start,
                              GridCell goal,
                              bool allow_unknown,
                              std::int8_t value = 0) {
  SmallGridCase property_case;
  property_case.map.header.frame_id = "map";
  property_case.map.info.width = width;
  property_case.map.info.height = height;
  property_case.map.info.resolution = 1.0F;
  property_case.map.info.origin.orientation.w = 1.0;
  property_case.map.data.assign(static_cast<std::size_t>(width) * height, value);
  property_case.start = start;
  property_case.goal = goal;
  property_case.allow_unknown = allow_unknown;
  return property_case;
}

std::size_t flatIndex(const SmallGridCase& property_case, GridCell cell) {
  return static_cast<std::size_t>(cell.y) * property_case.map.info.width + cell.x;
}

bool generatedCasesEqual(const SmallGridCase& lhs, const SmallGridCase& rhs) {
  return lhs.case_index == rhs.case_index && lhs.seed == rhs.seed && lhs.map == rhs.map &&
         lhs.start == rhs.start && lhs.goal == rhs.goal && lhs.allow_unknown == rhs.allow_unknown &&
         lhs.expected_reachable == rhs.expected_reachable;
}

TEST(SmallGridProperty, SeedScheduleIsStableAndUnique) {
  std::set<std::uint64_t> seeds;
  for (std::size_t case_index = 0; case_index < kCorpusSize; ++case_index) {
    const std::uint64_t seed = propertySeed(case_index);
    EXPECT_TRUE(seeds.insert(seed).second) << "case " << case_index;
    EXPECT_EQ(seed, kPropertyBaseSeed + kSeedStride * static_cast<std::uint64_t>(case_index + 1));
  }
}

TEST(SmallGridProperty, FixedSeedCorpusIsExactlyReproducible) {
  for (std::size_t case_index = 0; case_index < kCorpusSize; ++case_index) {
    const SmallGridCase first = makeSmallGridCase(case_index);
    SCOPED_TRACE(describeSmallGridCase(first));
    const SmallGridCase second = makeSmallGridCase(case_index);
    EXPECT_TRUE(generatedCasesEqual(first, second));
  }
}

TEST(SmallGridProperty, FixedSeedCorpusSatisfiesOraclePreconditions) {
  std::array<std::size_t, 4> combination_counts{};
  std::size_t unknown_values = 0;
  std::size_t known_occupied_values = 0;

  for (std::size_t case_index = 0; case_index < kCorpusSize; ++case_index) {
    const SmallGridCase property_case = makeSmallGridCase(case_index);
    SCOPED_TRACE(describeSmallGridCase(property_case));
    const auto& map = property_case.map;

    ASSERT_GE(map.info.width, 3U);
    ASSERT_GE(map.info.height, 3U);
    ASSERT_EQ(map.data.size(), static_cast<std::size_t>(map.info.width) * map.info.height);
    EXPECT_EQ(map.header.frame_id, "map");
    EXPECT_GT(map.info.resolution, 0.0F);
    EXPECT_EQ(map.info.origin.orientation.w, 1.0);

    ASSERT_LT(property_case.start.x, map.info.width);
    ASSERT_LT(property_case.start.y, map.info.height);
    ASSERT_LT(property_case.goal.x, map.info.width);
    ASSERT_LT(property_case.goal.y, map.info.height);
    EXPECT_FALSE(isForcedOuterRing(map, property_case.start.x, property_case.start.y));
    EXPECT_FALSE(isForcedOuterRing(map, property_case.goal.x, property_case.goal.y));
    EXPECT_FALSE(isEffectivelyBlocked(
      map, property_case.start.x, property_case.start.y, property_case.allow_unknown));
    EXPECT_FALSE(isEffectivelyBlocked(
      map, property_case.goal.x, property_case.goal.y, property_case.allow_unknown));
    EXPECT_FALSE(property_case.start == property_case.goal);

    for (const std::int8_t raw_value : map.data) {
      const int value = static_cast<int>(raw_value);
      EXPECT_TRUE(value == -1 || (value >= 0 && value <= 100))
        << "invalid raw OccupancyGrid value " << value;
      if (value == -1)
        ++unknown_values;
      else if (value >= kOccupiedThreshold)
        ++known_occupied_values;
    }

    EXPECT_FALSE(hasUnsupportedDiagonalContact(property_case));
    EXPECT_EQ(fourNeighborReachable(property_case), property_case.expected_reachable);

    const std::size_t combination =
      (property_case.allow_unknown ? 2U : 0U) + (property_case.expected_reachable ? 1U : 0U);
    ASSERT_LT(combination, combination_counts.size());
    ++combination_counts[combination];
  }

  EXPECT_EQ(combination_counts, (std::array<std::size_t, 4>{6, 6, 6, 6}));
  EXPECT_GT(unknown_values, 0U);
  EXPECT_GT(known_occupied_values, 0U);
}

TEST(SmallGridProperty, ForcedOuterRingIsBlockedEvenWhenRawCellsAreFree) {
  const SmallGridCase property_case = makeUniformCase(5, 5, GridCell{1, 1}, GridCell{3, 3}, false);
  ASSERT_TRUE(fourNeighborReachable(property_case));

  for (std::uint32_t y = 0; y < property_case.map.info.height; ++y) {
    for (std::uint32_t x = 0; x < property_case.map.info.width; ++x) {
      if (!isForcedOuterRing(property_case.map, x, y))
        continue;
      EXPECT_EQ(property_case.map.data[flatIndex(property_case, GridCell{x, y})], 0);
      EXPECT_TRUE(isEffectivelyBlocked(property_case.map, x, y, property_case.allow_unknown));
    }
  }

  SmallGridCase ring_endpoint = property_case;
  ring_endpoint.start = GridCell{0, 2};
  EXPECT_FALSE(fourNeighborReachable(ring_endpoint));
  ring_endpoint.start = property_case.start;
  ring_endpoint.goal = GridCell{4, 2};
  EXPECT_FALSE(fourNeighborReachable(ring_endpoint));
}

TEST(SmallGridProperty, OccupiedThresholdUsesGreaterThanOrEqualSemantics) {
  SmallGridCase property_case = makeUniformCase(7, 5, GridCell{1, 2}, GridCell{5, 2}, false);
  for (std::uint32_t y = 1; y + 1 < property_case.map.info.height; ++y)
    property_case.map.data[flatIndex(property_case, GridCell{3, y})] = 98;

  EXPECT_TRUE(fourNeighborReachable(property_case, 99));
  EXPECT_FALSE(fourNeighborReachable(property_case, 98));

  property_case.map.data[flatIndex(property_case, GridCell{3, 2})] = 99;
  EXPECT_TRUE(isEffectivelyBlocked(property_case.map, 3, 2, false, kOccupiedThreshold));
  EXPECT_FALSE(isEffectivelyBlocked(property_case.map, 3, 1, false, kOccupiedThreshold));

  property_case.map.data[flatIndex(property_case, GridCell{2, 2})] = -2;
  EXPECT_FALSE(hasValidOccupancyValues(property_case.map));
  EXPECT_FALSE(fourNeighborReachable(property_case));
  property_case.map.data[flatIndex(property_case, GridCell{2, 2})] = 101;
  EXPECT_FALSE(hasValidOccupancyValues(property_case.map));
  EXPECT_FALSE(fourNeighborReachable(property_case));
}

TEST(SmallGridProperty, UnknownCellsFollowAllowUnknownPolicy) {
  SmallGridCase property_case = makeUniformCase(7, 5, GridCell{1, 2}, GridCell{5, 2}, false);
  for (std::uint32_t y = 1; y + 1 < property_case.map.info.height; ++y)
    property_case.map.data[flatIndex(property_case, GridCell{3, y})] = -1;

  EXPECT_TRUE(isEffectivelyBlocked(property_case.map, 3, 2, false));
  EXPECT_FALSE(fourNeighborReachable(property_case));

  property_case.allow_unknown = true;
  EXPECT_FALSE(isEffectivelyBlocked(property_case.map, 3, 2, true));
  EXPECT_TRUE(fourNeighborReachable(property_case));
}

TEST(SmallGridProperty, DiagonalOnlyContactDetectorRecognizesSharedVertex) {
  SmallGridCase property_case = makeUniformCase(7, 7, GridCell{1, 1}, GridCell{5, 5}, false);
  property_case.map.data[flatIndex(property_case, GridCell{2, 2})] = 100;
  property_case.map.data[flatIndex(property_case, GridCell{3, 3})] = 100;

  EXPECT_TRUE(hasUnsupportedDiagonalContact(property_case));

  property_case.map.data[flatIndex(property_case, GridCell{2, 3})] = 100;
  EXPECT_FALSE(hasUnsupportedDiagonalContact(property_case))
    << "Three blocked cells in a 2x2 neighborhood are edge-connected";

  property_case = makeUniformCase(7, 7, GridCell{1, 1}, GridCell{5, 5}, false);
  property_case.map.data[flatIndex(property_case, GridCell{2, 2})] = 100;
  property_case.map.data[flatIndex(property_case, GridCell{3, 2})] = 100;
  EXPECT_FALSE(hasUnsupportedDiagonalContact(property_case))
    << "Edge-adjacent blocked cells do not create a shared-only vertex";

  for (std::uint32_t y = 0; y < property_case.map.info.height; ++y)
    property_case.map.data[flatIndex(property_case, GridCell{3, y})] = 100;
  EXPECT_FALSE(hasUnsupportedDiagonalContact(property_case))
    << "A barrier may join the forced outer ring along complete edges";
}

}  // namespace
