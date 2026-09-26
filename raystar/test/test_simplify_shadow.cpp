// Shadow-mode regression gate (review finding 4): with
// RAYSTAR_SIMPLIFY_SHADOW=1 every removability verdict is recomputed
// against a full-enumeration candidate set and any divergence aborts the
// process.  The flag is latched on the FIRST simplify call in the process,
// so this test lives in its own binary and sets the variable before any
// Polymap is built.

#include <gtest/gtest.h>
#include <raystar/polymap.h>

#include <cstdint>
#include <cstdlib>

#include "simplify_test_worlds.h"

namespace raystar {
namespace {

using simplify_test_worlds::connected;
using simplify_test_worlds::randomWorld;
using simplify_test_worlds::Rng;

TEST(SimplifyShadow, VerdictsAgreeOnRandomWorlds) {
  setenv("RAYSTAR_SIMPLIFY_SHADOW", "1", 1);
  int worlds = 0;
  for (std::uint64_t seed = 1; seed <= 8 && worlds < 5; ++seed) {
    Rng rng{20260926ULL * 17 + seed};
    const int width = rng.range(30, 42);
    const int height = rng.range(26, 36);
    const GridMap map = randomWorld(rng, width, height);
    int start_x = -1, start_y = -1, goal_x = -1, goal_y = -1;
    for (int attempt = 0; attempt < 60 && start_x < 0; ++attempt) {
      const int sx = rng.range(1, width - 2), sy = rng.range(1, height - 2);
      const int gx = rng.range(1, width - 2), gy = rng.range(1, height - 2);
      if ((sx == gx && sy == gy) ||
          map.data[static_cast<size_t>(sy) * map.width + static_cast<size_t>(sx)] != 0 ||
          map.data[static_cast<size_t>(gy) * map.width + static_cast<size_t>(gx)] != 0)
        continue;
      if (connected(map, sx, sy, gx, gy)) {
        start_x = sx;
        start_y = sy;
        goal_x = gx;
        goal_y = gy;
      }
    }
    if (start_x < 0)
      continue;
    const Point2d start{start_x + 0.5, start_y + 0.5};
    const Point2d goal{goal_x + 0.5, goal_y + 0.5};
    auto created = Polymap::create(
      map, start_x, start_y, start, {PolymapEndpoint{goal_x, goal_y, goal}}, StopToken{});
    if (created)
      ++worlds;
  }
  EXPECT_GE(worlds, 3);  // divergence would have aborted before this line
}

}  // namespace
}  // namespace raystar
