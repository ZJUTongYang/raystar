#pragma once

#include <raystar_interfaces/map_identity.hpp>

#include <nav_msgs/msg/occupancy_grid.hpp>
#include <unique_identifier_msgs/msg/uuid.hpp>

#include <cstdint>

namespace raystar_interfaces {

using EnvironmentId = unique_identifier_msgs::msg::UUID;

// These versions are part of the public environment-identity contract.  Bump
// the corresponding value whenever the interpretation changes, even when the
// OccupancyGrid bytes and request policy remain identical.
inline constexpr std::uint32_t kEnvironmentIdentityVersion = 1U;
inline constexpr std::uint32_t kOccupancySemanticsVersion = 1U;
inline constexpr std::uint32_t kGeometrySemanticsVersion = 1U;
inline constexpr std::uint32_t kTopologySemanticsVersion = 1U;

struct EnvironmentSemanticVersions {
  std::uint32_t identity{kEnvironmentIdentityVersion};
  std::uint32_t occupancy{kOccupancySemanticsVersion};
  std::uint32_t geometry{kGeometrySemanticsVersion};
  std::uint32_t topology{kTopologySemanticsVersion};
};

/// Compute a deterministic identity for the planning environment shared by
/// Raystar GCP and UPS.  Unlike map_id, this excludes message and map-load
/// timestamps: republishing identical planning content preserves the identity.
/// It includes every OccupancyGrid field interpreted by the 2-D planner, the
/// raw occupancy payload, the fixed threshold and per-request unknown policy,
/// and explicit semantic versions.  It is an identity guard, not a
/// cryptographic authenticity proof.
inline EnvironmentId computeEnvironmentId(
  const nav_msgs::msg::OccupancyGrid& map,
  std::int32_t occupied_threshold,
  bool allow_unknown,
  const EnvironmentSemanticVersions& versions = {}) noexcept {
  detail::MapIdentityHasher hasher;
  hasher.addString("raystar-planning-environment");
  hasher.addInteger(versions.identity);
  hasher.addInteger(versions.occupancy);
  hasher.addInteger(versions.geometry);
  hasher.addInteger(versions.topology);

  hasher.addString(map.header.frame_id);
  hasher.addFloating(map.info.resolution);
  hasher.addInteger(map.info.width);
  hasher.addInteger(map.info.height);
  hasher.addFloating(map.info.origin.position.x);
  hasher.addFloating(map.info.origin.position.y);
  hasher.addFloating(map.info.origin.position.z);
  hasher.addFloating(map.info.origin.orientation.x);
  hasher.addFloating(map.info.origin.orientation.y);
  hasher.addFloating(map.info.origin.orientation.z);
  hasher.addFloating(map.info.origin.orientation.w);
  hasher.addInteger<std::uint64_t>(static_cast<std::uint64_t>(map.data.size()));
  for (const std::int8_t value : map.data) hasher.addByte(static_cast<std::uint8_t>(value));

  hasher.addInteger(occupied_threshold);
  hasher.addByte(allow_unknown ? 1U : 0U);
  return hasher.finish();
}

inline bool environmentIdsEqual(const EnvironmentId& first, const EnvironmentId& second) noexcept {
  return first.uuid == second.uuid;
}

inline bool isZeroEnvironmentId(const EnvironmentId& id) noexcept {
  return isZeroMapId(id);
}

}  // namespace raystar_interfaces
