#ifndef RAYSTAR_TEST_SMALL_GRID_PROPERTY_H_
#define RAYSTAR_TEST_SMALL_GRID_PROPERTY_H_

#include <nav_msgs/msg/occupancy_grid.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// Test-only fixed-seed OccupancyGrid generation and a four-neighbor BFS
// reachability oracle.  This file deliberately does not include any Raystar
// production header, so the differential check cannot accidentally reuse the
// planner's map conversion, contour extraction, or reachability logic.
namespace raystar::test_property
{

inline constexpr int kOccupiedThreshold = 99;
inline constexpr std::uint64_t kPropertyBaseSeed = 0x5241595354415234ULL;
inline constexpr std::uint64_t kSeedStride = 0x9e3779b97f4a7c15ULL;

struct GridCell
{
  std::uint32_t x = 0;
  std::uint32_t y = 0;

  bool operator==(const GridCell & other) const noexcept
  {
    return x == other.x && y == other.y;
  }
};

struct SmallGridCase
{
  std::size_t case_index = 0;
  std::uint64_t seed = 0;
  nav_msgs::msg::OccupancyGrid map;
  GridCell start;
  GridCell goal;
  bool allow_unknown = false;
  bool expected_reachable = false;
};

class DeterministicRng
{
public:
  explicit DeterministicRng(std::uint64_t seed) : state_(seed) {}

  std::uint64_t next() noexcept
  {
    // SplitMix64 is specified here instead of using a standard-library
    // distribution, whose integer mapping is not required to be identical
    // across implementations.
    std::uint64_t value = (state_ += kSeedStride);
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
  }

  std::size_t bounded(std::size_t upper_exclusive) noexcept
  {
    return upper_exclusive == 0 ? 0 :
      static_cast<std::size_t>(next() % upper_exclusive);
  }

private:
  std::uint64_t state_;
};

inline std::uint64_t propertySeed(std::size_t case_index) noexcept
{
  return kPropertyBaseSeed +
    kSeedStride * static_cast<std::uint64_t>(case_index + 1);
}

inline bool isForcedOuterRing(
  const nav_msgs::msg::OccupancyGrid & map,
  std::uint32_t x, std::uint32_t y) noexcept
{
  return x == 0 || y == 0 || x + 1 == map.info.width ||
    y + 1 == map.info.height;
}

inline bool hasValidOccupancyValues(
  const nav_msgs::msg::OccupancyGrid & map) noexcept
{
  return std::all_of(
    map.data.begin(), map.data.end(), [](std::int8_t raw_value) {
      const int value = static_cast<int>(raw_value);
      return value == -1 || (value >= 0 && value <= 100);
    });
}

inline bool isEffectivelyBlocked(
  const nav_msgs::msg::OccupancyGrid & map,
  std::uint32_t x, std::uint32_t y, bool allow_unknown,
  int occupied_threshold = kOccupiedThreshold) noexcept
{
  if (x >= map.info.width || y >= map.info.height ||
      isForcedOuterRing(map, x, y))
  {
    return true;
  }
  const std::size_t index =
    static_cast<std::size_t>(y) * map.info.width + x;
  if (index >= map.data.size())
    return true;
  const int value = static_cast<int>(map.data[index]);
  if (value == -1)
    return !allow_unknown;
  return value >= occupied_threshold;
}

inline bool fourNeighborReachable(
  const SmallGridCase & property_case,
  int occupied_threshold = kOccupiedThreshold)
{
  const auto & map = property_case.map;
  if (map.info.width == 0 || map.info.height == 0 ||
      map.data.size() !=
      static_cast<std::size_t>(map.info.width) * map.info.height ||
      property_case.start.x >= map.info.width ||
      property_case.start.y >= map.info.height ||
      property_case.goal.x >= map.info.width ||
      property_case.goal.y >= map.info.height ||
      occupied_threshold < 1 || occupied_threshold > 100 ||
      !hasValidOccupancyValues(map))
  {
    return false;
  }
  if (isEffectivelyBlocked(
      map, property_case.start.x, property_case.start.y,
      property_case.allow_unknown, occupied_threshold) ||
      isEffectivelyBlocked(
      map, property_case.goal.x, property_case.goal.y,
      property_case.allow_unknown, occupied_threshold))
  {
    return false;
  }

  const std::size_t cell_count =
    static_cast<std::size_t>(map.info.width) * map.info.height;
  std::vector<std::uint8_t> visited(cell_count, 0);
  std::vector<GridCell> frontier;
  frontier.reserve(cell_count);
  const auto flat_index = [&](const GridCell & cell) {
    return static_cast<std::size_t>(cell.y) * map.info.width + cell.x;
  };
  frontier.push_back(property_case.start);
  visited[flat_index(property_case.start)] = 1;

  constexpr std::array<int, 4> dx = {1, -1, 0, 0};
  constexpr std::array<int, 4> dy = {0, 0, 1, -1};
  for (std::size_t head = 0; head < frontier.size(); ++head)
  {
    const GridCell current = frontier[head];
    if (current == property_case.goal)
      return true;
    for (std::size_t direction = 0; direction < dx.size(); ++direction)
    {
      const int next_x = static_cast<int>(current.x) + dx[direction];
      const int next_y = static_cast<int>(current.y) + dy[direction];
      if (next_x < 0 || next_y < 0 ||
          next_x >= static_cast<int>(map.info.width) ||
          next_y >= static_cast<int>(map.info.height))
      {
        continue;
      }
      const GridCell next{
        static_cast<std::uint32_t>(next_x),
        static_cast<std::uint32_t>(next_y)};
      const std::size_t next_index = flat_index(next);
      if (visited[next_index] || isEffectivelyBlocked(
          map, next.x, next.y, property_case.allow_unknown,
          occupied_threshold))
      {
        continue;
      }
      visited[next_index] = 1;
      frontier.push_back(next);
    }
  }
  return false;
}

// The contour extractor deliberately rejects two blocked components that meet
// only at one grid vertex.  T-junction/crossing contour relations cannot arise
// from the generator's separated rectangles once this diagonal pattern is
// excluded.  A full barrier is allowed to join the forced outer ring by edges.
inline bool hasUnsupportedDiagonalContact(
  const SmallGridCase & property_case,
  int occupied_threshold = kOccupiedThreshold)
{
  const auto & map = property_case.map;
  if (map.info.width < 2 || map.info.height < 2)
    return false;
  for (std::uint32_t y = 0; y + 1 < map.info.height; ++y)
  {
    for (std::uint32_t x = 0; x + 1 < map.info.width; ++x)
    {
      const bool lower_left = isEffectivelyBlocked(
        map, x, y, property_case.allow_unknown, occupied_threshold);
      const bool lower_right = isEffectivelyBlocked(
        map, x + 1, y, property_case.allow_unknown, occupied_threshold);
      const bool upper_left = isEffectivelyBlocked(
        map, x, y + 1, property_case.allow_unknown, occupied_threshold);
      const bool upper_right = isEffectivelyBlocked(
        map, x + 1, y + 1, property_case.allow_unknown,
        occupied_threshold);
      if ((lower_left && upper_right && !lower_right && !upper_left) ||
          (lower_right && upper_left && !lower_left && !upper_right))
      {
        return true;
      }
    }
  }
  return false;
}

inline std::int8_t randomKnownFreeValue(DeterministicRng & rng)
{
  constexpr std::array<std::int8_t, 5> values = {0, 1, 25, 50, 98};
  return values[rng.bounded(values.size())];
}

inline SmallGridCase makeSmallGridCase(std::size_t case_index)
{
  SmallGridCase property_case;
  property_case.case_index = case_index;
  property_case.seed = propertySeed(case_index);
  // Cycle through all four (reachable/unreachable) x (unknown blocked/free)
  // combinations while keeping a balanced 12/12 differential corpus.
  property_case.allow_unknown = (case_index / 2) % 2 == 1;
  property_case.expected_reachable = case_index % 2 == 1;

  DeterministicRng rng(property_case.seed);
  auto & map = property_case.map;
  map.header.frame_id = "map";
  map.info.width = 9 + static_cast<std::uint32_t>(rng.bounded(5));
  map.info.height = 9 + static_cast<std::uint32_t>(rng.bounded(5));
  constexpr std::array<float, 4> resolutions = {0.25F, 0.5F, 1.0F, 2.0F};
  map.info.resolution = resolutions[rng.bounded(resolutions.size())];
  map.info.origin.position.x =
    0.25 * (static_cast<int>(rng.bounded(17)) - 8);
  map.info.origin.position.y =
    0.25 * (static_cast<int>(rng.bounded(17)) - 8);
  map.info.origin.orientation.w = 1.0;

  const std::size_t cell_count =
    static_cast<std::size_t>(map.info.width) * map.info.height;
  std::vector<std::uint8_t> blocked(cell_count, 0);
  const auto mask_index = [&](std::uint32_t x, std::uint32_t y) {
    return static_cast<std::size_t>(y) * map.info.width + x;
  };

  if (!property_case.expected_reachable)
  {
    // A one-cell barrier joins opposite sides of the forced outer ring.  No
    // other obstacle is added, keeping the contour topology supported while
    // producing two unambiguously disconnected free-space components.
    const bool vertical = rng.bounded(2) == 0;
    if (vertical)
    {
      const std::uint32_t barrier_x = 3 +
        static_cast<std::uint32_t>(rng.bounded(map.info.width - 6));
      for (std::uint32_t y = 0; y < map.info.height; ++y)
        blocked[mask_index(barrier_x, y)] = 1;
      property_case.start = GridCell{
        1 + static_cast<std::uint32_t>(rng.bounded(barrier_x - 1)),
        1 + static_cast<std::uint32_t>(rng.bounded(map.info.height - 2))};
      property_case.goal = GridCell{
        barrier_x + 1 + static_cast<std::uint32_t>(
          rng.bounded(map.info.width - barrier_x - 2)),
        1 + static_cast<std::uint32_t>(rng.bounded(map.info.height - 2))};
    }
    else
    {
      const std::uint32_t barrier_y = 3 +
        static_cast<std::uint32_t>(rng.bounded(map.info.height - 6));
      for (std::uint32_t x = 0; x < map.info.width; ++x)
        blocked[mask_index(x, barrier_y)] = 1;
      property_case.start = GridCell{
        1 + static_cast<std::uint32_t>(rng.bounded(map.info.width - 2)),
        1 + static_cast<std::uint32_t>(rng.bounded(barrier_y - 1))};
      property_case.goal = GridCell{
        1 + static_cast<std::uint32_t>(rng.bounded(map.info.width - 2)),
        barrier_y + 1 + static_cast<std::uint32_t>(
          rng.bounded(map.info.height - barrier_y - 2))};
    }
  }
  else
  {
    // Place up to four 1x1..3x3 rectangles.  Every candidate stays one free
    // cell away from the forced ring and the other rectangles, excluding both
    // shared edges and diagonal-only shared contour vertices.
    const std::size_t target_rectangles = rng.bounded(5);
    std::size_t accepted_rectangles = 0;
    for (std::size_t attempt = 0;
         attempt < 64 && accepted_rectangles < target_rectangles; ++attempt)
    {
      const std::uint32_t left = 2 +
        static_cast<std::uint32_t>(rng.bounded(map.info.width - 4));
      const std::uint32_t bottom = 2 +
        static_cast<std::uint32_t>(rng.bounded(map.info.height - 4));
      const std::uint32_t right = std::min<std::uint32_t>(
        map.info.width - 3,
        left + static_cast<std::uint32_t>(rng.bounded(3)));
      const std::uint32_t top = std::min<std::uint32_t>(
        map.info.height - 3,
        bottom + static_cast<std::uint32_t>(rng.bounded(3)));

      bool separated = true;
      for (std::uint32_t y = bottom - 1; y <= top + 1 && separated; ++y)
      {
        for (std::uint32_t x = left - 1; x <= right + 1; ++x)
        {
          if (blocked[mask_index(x, y)])
          {
            separated = false;
            break;
          }
        }
      }
      if (!separated)
        continue;
      for (std::uint32_t y = bottom; y <= top; ++y)
        for (std::uint32_t x = left; x <= right; ++x)
          blocked[mask_index(x, y)] = 1;
      ++accepted_rectangles;
    }

    std::vector<GridCell> free_cells;
    for (std::uint32_t y = 1; y + 1 < map.info.height; ++y)
    {
      for (std::uint32_t x = 1; x + 1 < map.info.width; ++x)
      {
        if (!blocked[mask_index(x, y)])
          free_cells.push_back(GridCell{x, y});
      }
    }
    // The dimensions and one-cell rectangle moats above guarantee many free
    // cells.  Keep a defensive fallback here nevertheless, so future changes
    // to the generator cannot turn the index arithmetic below into an
    // underflow or out-of-bounds access.
    if (free_cells.size() < 2)
    {
      std::fill(blocked.begin(), blocked.end(), 0);
      free_cells = {
        GridCell{1, 1},
        GridCell{map.info.width - 2, map.info.height - 2}};
    }
    const std::size_t start_index = rng.bounded(free_cells.size());
    const std::size_t goal_index =
      (start_index + 1 + rng.bounded(free_cells.size() - 1)) %
      free_cells.size();
    property_case.start = free_cells[start_index];
    property_case.goal = free_cells[goal_index];
  }

  map.data.resize(cell_count);
  for (std::uint32_t y = 0; y < map.info.height; ++y)
  {
    for (std::uint32_t x = 0; x < map.info.width; ++x)
    {
      const std::size_t index = mask_index(x, y);
      if (blocked[index])
      {
        map.data[index] = !property_case.allow_unknown && rng.bounded(4) == 0 ?
          static_cast<std::int8_t>(-1) :
          static_cast<std::int8_t>(99 + rng.bounded(2));
      }
      else
      {
        map.data[index] = property_case.allow_unknown && rng.bounded(5) == 0 ?
          static_cast<std::int8_t>(-1) : randomKnownFreeValue(rng);
      }
    }
  }
  return property_case;
}

inline std::string describeSmallGridCase(const SmallGridCase & property_case)
{
  const auto & map = property_case.map;
  std::ostringstream stream;
  stream << "T-04 case=" << property_case.case_index << " seed=0x"
         << std::hex << std::uppercase << property_case.seed << std::dec
         << " size=" << map.info.width << 'x' << map.info.height
         << " resolution=" << map.info.resolution
         << " origin=(" << map.info.origin.position.x << ','
         << map.info.origin.position.y << ')'
         << " allow_unknown=" << std::boolalpha
         << property_case.allow_unknown
         << " expected_reachable=" << property_case.expected_reachable
         << " start=(" << property_case.start.x << ','
         << property_case.start.y << ") goal=(" << property_case.goal.x
         << ',' << property_case.goal.y << ")\n"
         << "effective grid (top row first; X=forced ring, #=blocked, "
            "?=unknown, .=known free):\n";

  for (std::uint32_t row = map.info.height; row > 0; --row)
  {
    const std::uint32_t y = row - 1;
    stream << std::setw(2) << y << ' ';
    for (std::uint32_t x = 0; x < map.info.width; ++x)
    {
      char symbol = '.';
      if (GridCell{x, y} == property_case.start)
        symbol = 'S';
      else if (GridCell{x, y} == property_case.goal)
        symbol = 'G';
      else if (isForcedOuterRing(map, x, y))
        symbol = 'X';
      else
      {
        const int value = static_cast<int>(
          map.data[static_cast<std::size_t>(y) * map.info.width + x]);
        if (value == -1)
          symbol = '?';
        else if (value >= kOccupiedThreshold)
          symbol = '#';
      }
      stream << symbol;
    }
    stream << '\n';
  }

  stream << "raw OccupancyGrid rows (top row first):\n";
  for (std::uint32_t row = map.info.height; row > 0; --row)
  {
    const std::uint32_t y = row - 1;
    stream << std::setw(2) << y << ':';
    for (std::uint32_t x = 0; x < map.info.width; ++x)
    {
      stream << ' ' << std::setw(3) << static_cast<int>(
        map.data[static_cast<std::size_t>(y) * map.info.width + x]);
    }
    stream << '\n';
  }
  return stream.str();
}

}  // namespace raystar::test_property

#endif  // RAYSTAR_TEST_SMALL_GRID_PROPERTY_H_
