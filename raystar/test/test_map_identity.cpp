#include <gtest/gtest.h>

#include <raystar_interfaces/environment_identity.hpp>
#include <raystar_interfaces/map_identity.hpp>

#include <algorithm>
#include <cstdint>

namespace {

nav_msgs::msg::OccupancyGrid makeMap() {
  nav_msgs::msg::OccupancyGrid map;
  map.header.frame_id = "map";
  map.header.stamp.sec = 12;
  map.header.stamp.nanosec = 34;
  map.info.map_load_time.sec = 5;
  map.info.map_load_time.nanosec = 6;
  map.info.resolution = 0.25F;
  map.info.width = 3;
  map.info.height = 2;
  map.info.origin.position.x = -1.0;
  map.info.origin.position.y = 2.0;
  map.info.origin.orientation.w = 1.0;
  map.data = {0, 10, 100, -1, 0, 50};
  return map;
}

TEST(MapIdentity, IsDeterministicAndUsesApplicationUuidBits) {
  const auto map = makeMap();
  const auto first = raystar_interfaces::computeMapId(map);
  const auto second = raystar_interfaces::computeMapId(map);

  EXPECT_TRUE(raystar_interfaces::mapIdsEqual(first, second));
  EXPECT_FALSE(raystar_interfaces::isZeroMapId(first));
  EXPECT_EQ(first.uuid[6] >> 4U, 8U);
  EXPECT_EQ(first.uuid[8] & 0xc0U, 0x80U);
}

TEST(MapIdentity, EveryPlanningRelevantFieldChangesIdentity) {
  const auto baseline_map = makeMap();
  const auto baseline = raystar_interfaces::computeMapId(baseline_map);

  const auto expect_changed = [&baseline](const auto& candidate) {
    EXPECT_FALSE(
      raystar_interfaces::mapIdsEqual(baseline, raystar_interfaces::computeMapId(candidate)));
  };

  auto candidate = baseline_map;
  candidate.header.frame_id = "other";
  expect_changed(candidate);
  candidate = baseline_map;
  ++candidate.header.stamp.nanosec;
  expect_changed(candidate);
  candidate = baseline_map;
  ++candidate.info.map_load_time.sec;
  expect_changed(candidate);
  candidate = baseline_map;
  candidate.info.resolution = 0.5F;
  expect_changed(candidate);
  candidate = baseline_map;
  ++candidate.info.width;
  expect_changed(candidate);
  candidate = baseline_map;
  candidate.info.origin.position.x += 1.0;
  expect_changed(candidate);
  candidate = baseline_map;
  candidate.info.origin.orientation.w = -1.0;
  expect_changed(candidate);
  candidate = baseline_map;
  candidate.data.back() = 51;
  expect_changed(candidate);
  candidate = baseline_map;
  candidate.data.push_back(0);
  expect_changed(candidate);
}

TEST(MapIdentity, ZeroIdentityIsReservedForMissingReferences) {
  raystar_interfaces::MapId zero;
  EXPECT_TRUE(raystar_interfaces::isZeroMapId(zero));

  auto non_zero = zero;
  non_zero.uuid.back() = 1U;
  EXPECT_FALSE(raystar_interfaces::isZeroMapId(non_zero));
  EXPECT_FALSE(raystar_interfaces::mapIdsEqual(zero, non_zero));
}

TEST(EnvironmentIdentity, IsStableAcrossTimestampOnlyRepublishes) {
  const auto baseline_map = makeMap();
  auto republished_map = baseline_map;
  ++republished_map.header.stamp.sec;
  ++republished_map.info.map_load_time.nanosec;

  const auto baseline =
    raystar_interfaces::computeEnvironmentId(baseline_map, 99, false);
  const auto republished =
    raystar_interfaces::computeEnvironmentId(republished_map, 99, false);

  EXPECT_TRUE(raystar_interfaces::environmentIdsEqual(baseline, republished));
  EXPECT_FALSE(raystar_interfaces::isZeroEnvironmentId(baseline));
  EXPECT_FALSE(raystar_interfaces::mapIdsEqual(
    raystar_interfaces::computeMapId(baseline_map),
    raystar_interfaces::computeMapId(republished_map)));
}

TEST(EnvironmentIdentity, BindsMapPolicyAndSemanticVersions) {
  const auto map = makeMap();
  const auto baseline = raystar_interfaces::computeEnvironmentId(map, 99, false);
  const auto expect_changed = [&baseline](const auto& candidate) {
    EXPECT_FALSE(raystar_interfaces::environmentIdsEqual(baseline, candidate));
  };

  auto changed_map = map;
  changed_map.data.back() = 51;
  expect_changed(raystar_interfaces::computeEnvironmentId(changed_map, 99, false));
  changed_map = map;
  changed_map.info.origin.position.x += 0.5;
  expect_changed(raystar_interfaces::computeEnvironmentId(changed_map, 99, false));
  expect_changed(raystar_interfaces::computeEnvironmentId(map, 50, false));
  expect_changed(raystar_interfaces::computeEnvironmentId(map, 99, true));

  raystar_interfaces::EnvironmentSemanticVersions versions;
  ++versions.occupancy;
  expect_changed(raystar_interfaces::computeEnvironmentId(map, 99, false, versions));
  versions = {};
  ++versions.geometry;
  expect_changed(raystar_interfaces::computeEnvironmentId(map, 99, false, versions));
  versions = {};
  ++versions.topology;
  expect_changed(raystar_interfaces::computeEnvironmentId(map, 99, false, versions));
}

}  // namespace
