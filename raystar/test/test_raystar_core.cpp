#include <gtest/gtest.h>
#include <raystar/raystar_core.h>

using namespace raystar;

static GridMap makeCorridorMap()
{
  GridMap map;
  map.width = 20;
  map.height = 30;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.resize(600, 0);

  for (unsigned int x = 8; x < 12; ++x)
    for (unsigned int y = 0; y < 30; ++y)
      map.data[y * 20 + x] = (y < 8 || y > 11) ? 1 : 0;

  return map;
}

static GridMap makeSimpleMap()
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

TEST(RaystarCore, PlanStraightPath)
{
  GridMap map;
  map.width = 20;
  map.height = 20;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.resize(400, 0);

  RaystarCore core;
  auto result = core.plan(map, 2, 2, 17, 17, 1, false);

  EXPECT_TRUE(result.success);
  EXPECT_GE(result.path_solutions.size(), 1u);
  EXPECT_GT(result.path_solutions[0].path_cost_, 0.0);
}

TEST(RaystarCore, PlanAroundObstacle)
{
  auto map = makeSimpleMap();
  RaystarCore core;
  auto result = core.plan(map, 5, 15, 25, 15, 1, false);

  EXPECT_TRUE(result.success);
  EXPECT_GE(result.path_solutions.size(), 1u);
}

TEST(RaystarCore, NoPathToGoal)
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

  RaystarCore core;
  auto result = core.plan(map, 5, 10, 15, 10, 1, false);

  EXPECT_FALSE(result.success);
}

TEST(RaystarCore, MultiplePaths)
{
  auto map = makeSimpleMap();
  RaystarCore core;
  auto result = core.plan(map, 5, 15, 25, 15, 5, false);

  EXPECT_TRUE(result.success);
  EXPECT_GE(result.path_solutions.size(), 1u);
}

TEST(RaystarCore, SelfCrossingPruning)
{
  auto map = makeSimpleMap();
  RaystarCore core;

  auto result_with = core.plan(map, 5, 15, 25, 15, 3, true);
  auto result_without = core.plan(map, 5, 15, 25, 15, 3, false);

  EXPECT_TRUE(result_with.success);
  EXPECT_TRUE(result_without.success);
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
