#pragma once

// Shared randomized-world helpers for the simplification property tests
// (test_simplify_coverage.cpp, test_simplify_shadow.cpp).

#include <raystar/polymap.h>

#include <cstdint>
#include <vector>

namespace raystar {
namespace simplify_test_worlds {

struct Rng {
  std::uint64_t state;
  std::uint32_t next() {
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;
    return static_cast<std::uint32_t>(state >> 33);
  }
  int range(int lo, int hi) {  // inclusive
    return lo + static_cast<int>(next() % static_cast<std::uint32_t>(hi - lo + 1));
  }
};

GridMap randomWorld(Rng& rng, int width, int height) {
  GridMap map;
  map.width = static_cast<unsigned int>(width);
  map.height = static_cast<unsigned int>(height);
  map.resolution = 1.0F;
  map.data.assign(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
  const auto block = [&](int x0, int y0, int x1, int y1) {
    for (int y = y0; y <= y1; ++y)
      for (int x = x0; x <= x1; ++x)
        if (x >= 0 && y >= 0 && x < width && y < height)
          map.data[static_cast<size_t>(y) * static_cast<size_t>(width) + static_cast<size_t>(x)] =
            1;
  };
  block(0, 0, width - 1, 0);
  block(0, height - 1, width - 1, height - 1);
  block(0, 0, 0, height - 1);
  block(width - 1, 0, width - 1, height - 1);
  const int walls = rng.range(2, 5);
  for (int wall = 0; wall < walls; ++wall) {
    if (rng.range(0, 1) == 1) {
      const int y = rng.range(4, height - 5);
      const int x0 = rng.range(2, width - 12);
      block(x0, y, x0 + rng.range(6, 14), y);
    } else {
      const int x = rng.range(4, width - 5);
      const int y0 = rng.range(2, height - 12);
      block(x, y0, x, y0 + rng.range(6, 14));
    }
  }
  const int blobs = rng.range(1, 4);
  for (int blob = 0; blob < blobs; ++blob) {
    const int x = rng.range(2, width - 6);
    const int y = rng.range(2, height - 6);
    block(x, y, x + rng.range(1, 3), y + rng.range(1, 3));
  }
  return map;
}

bool connected(const GridMap& map, int sx, int sy, int gx, int gy) {
  const int width = static_cast<int>(map.width);
  const int height = static_cast<int>(map.height);
  std::vector<char> seen(map.data.size(), 0);
  std::vector<int> stack{sy * width + sx};
  seen[static_cast<size_t>(sy * width + sx)] = 1;
  while (!stack.empty()) {
    const int cell = stack.back();
    stack.pop_back();
    if (cell == gy * width + gx)
      return true;
    const int x = cell % width;
    const int y = cell / width;
    const int neighbours[4][2] = {{x - 1, y}, {x + 1, y}, {x, y - 1}, {x, y + 1}};
    for (const auto& neighbour : neighbours) {
      if (neighbour[0] < 0 || neighbour[1] < 0 || neighbour[0] >= width || neighbour[1] >= height)
        continue;
      const int index = neighbour[1] * width + neighbour[0];
      if (map.data[static_cast<size_t>(index)] != 0 || seen[static_cast<size_t>(index)])
        continue;
      seen[static_cast<size_t>(index)] = 1;
      stack.push_back(index);
    }
  }
  return false;
}

}  // namespace simplify_test_worlds
}  // namespace raystar
