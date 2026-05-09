#include <gtest/gtest.h>
#include <raystar/polymap.h>

using namespace raystar;

static GridMap makeMapWithObstacle()
{
  GridMap map;
  map.width = 30;
  map.height = 30;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.resize(900, 0);

  for (unsigned int y = 10; y < 20; ++y)
    for (unsigned int x = 10; x < 20; ++x)
      map.data[y * 30 + x] = 1;

  return map;
}

TEST(Polymap, ConstructionDetectsObstacles)
{
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);

  EXPECT_TRUE(poly.solution_exist_);
  EXPECT_GE(poly.obs_.size(), 1u);
}

TEST(Polymap, NoSolutionWhenGoalBlocked)
{
  GridMap map;
  map.width = 20;
  map.height = 20;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.resize(400, 0);

  for (unsigned int y = 0; y < 20; ++y)
    map.data[y * 20 + 10] = 1;

  Polymap poly(map, 5, 10, 15, 10);
  EXPECT_FALSE(poly.solution_exist_);
}

TEST(Polymap, VerticesRegistered)
{
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);

  if (!poly.obs_.empty() && !poly.obs_[0].detailed_ordered_vertices_.empty()) {
    auto first = poly.obs_[0].detailed_ordered_vertices_[0];
    auto loc = poly.locateVertex(first.first, first.second);
    EXPECT_EQ(loc.first, 0);
    EXPECT_EQ(loc.second, 0);
  }
}

TEST(Polymap, VisibilityRegionComputed)
{
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);

  if (poly.solution_exist_) {
    EXPECT_GT(poly.obs_.size(), 0u);
    EXPECT_GT(poly.obs_[0].detailed_ordered_vertices_.size(), 0u);
  }
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
