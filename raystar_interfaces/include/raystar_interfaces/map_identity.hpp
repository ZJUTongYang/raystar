#pragma once

#include <nav_msgs/msg/occupancy_grid.hpp>
#include <unique_identifier_msgs/msg/uuid.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <string>
#include <type_traits>

namespace raystar_interfaces
{

using MapId = unique_identifier_msgs::msg::UUID;

namespace detail
{

class MapIdentityHasher
{
public:
  void addByte(std::uint8_t value) noexcept
  {
    first_ ^= value;
    first_ *= 1099511628211ULL;

    second_ ^= static_cast<std::uint64_t>(value) + 0x9dULL;
    second_ *= 14029467366897019727ULL;
    second_ ^= second_ >> 29;
  }

  template<typename IntegerT>
  void addInteger(IntegerT value) noexcept
  {
    using UnsignedT = std::make_unsigned_t<IntegerT>;
    const UnsignedT unsigned_value = static_cast<UnsignedT>(value);
    for (std::size_t index = 0; index < sizeof(UnsignedT); ++index)
    {
      addByte(static_cast<std::uint8_t>(
        unsigned_value >> static_cast<unsigned int>(index * 8U)));
    }
  }

  template<typename FloatingT>
  void addFloating(FloatingT value) noexcept
  {
    static_assert(
      std::is_same<FloatingT, float>::value ||
      std::is_same<FloatingT, double>::value,
      "Only ROS floating-point field types are supported");
    using BitsT = std::conditional_t<
      std::is_same<FloatingT, float>::value, std::uint32_t, std::uint64_t>;
    BitsT bits = 0;
    std::memcpy(&bits, &value, sizeof(value));
    addInteger(bits);
  }

  void addString(const std::string& value) noexcept
  {
    addInteger<std::uint64_t>(static_cast<std::uint64_t>(value.size()));
    for (const unsigned char character : value)
      addByte(character);
  }

  [[nodiscard]] MapId finish() const noexcept
  {
    MapId result;
    for (std::size_t index = 0; index < 8; ++index)
    {
      const unsigned int shift = static_cast<unsigned int>((7U - index) * 8U);
      result.uuid[index] = static_cast<std::uint8_t>(first_ >> shift);
      result.uuid[index + 8] = static_cast<std::uint8_t>(second_ >> shift);
    }

    // RFC 9562 version 8 / variant 1 marks this as an application-defined
    // deterministic UUID. It is an identity guard, not a cryptographic
    // authenticity proof.
    result.uuid[6] = static_cast<std::uint8_t>((result.uuid[6] & 0x0fU) | 0x80U);
    result.uuid[8] = static_cast<std::uint8_t>((result.uuid[8] & 0x3fU) | 0x80U);
    return result;
  }

private:
  std::uint64_t first_{14695981039346656037ULL};
  std::uint64_t second_{0x84222325cbf29ce4ULL};
};

}  // namespace detail

/// Compute a deterministic identity over every OccupancyGrid field that can
/// affect planning or map provenance. The function performs no allocation and
/// is shared by the server, RViz panel, and tests so requests never need to
/// carry the full grid merely to identify a cached snapshot.
inline MapId computeMapId(const nav_msgs::msg::OccupancyGrid& map) noexcept
{
  detail::MapIdentityHasher hasher;
  hasher.addString("raystar-occupancy-grid-v1");
  hasher.addString(map.header.frame_id);
  hasher.addInteger(map.header.stamp.sec);
  hasher.addInteger(map.header.stamp.nanosec);
  hasher.addInteger(map.info.map_load_time.sec);
  hasher.addInteger(map.info.map_load_time.nanosec);
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
  for (const std::int8_t value : map.data)
    hasher.addByte(static_cast<std::uint8_t>(value));
  return hasher.finish();
}

inline bool mapIdsEqual(const MapId& first, const MapId& second) noexcept
{
  return first.uuid == second.uuid;
}

inline bool isZeroMapId(const MapId& id) noexcept
{
  for (const std::uint8_t byte : id.uuid)
  {
    if (byte != 0U)
      return false;
  }
  return true;
}

}  // namespace raystar_interfaces
