#include <gtest/gtest.h>
#include <raystar/coordinate_utils.h>

using namespace raystar;

TEST(CoordinateUtils, WorldToMapBasic)
{
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

TEST(CoordinateUtils, WorldToMapOutOfRange)
{
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

TEST(CoordinateUtils, MapToWorldAndBack)
{
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
      mapToWorld(map, mx, my, wx, wy);
      unsigned int mx2, my2;
      ASSERT_TRUE(worldToMap(map, wx, wy, mx2, my2)) << "mx=" << mx << " my=" << my;
      EXPECT_EQ(mx, mx2) << "wx=" << wx;
      EXPECT_EQ(my, my2) << "wy=" << wy;
    }
  }
}

TEST(CoordinateUtils, GridMapAccess)
{
  GridMap map;
  map.width = 10;
  map.height = 10;
  map.data.resize(100, 0);

  map(5, 5) = 1;
  EXPECT_EQ(map(5, 5), 1);
  EXPECT_EQ(map(0, 0), 0);
  EXPECT_EQ(map(9u, 9u), 0);
}

TEST(CoordinateUtils, GridMapOutOfBoundsReturnsOccupied)
{
  GridMap map;
  map.width = 10;
  map.height = 10;
  map.data.resize(100, 0);

  EXPECT_EQ(map.at(10u, 0u), 1u);
  EXPECT_EQ(map.at(0u, 10u), 1u);
}

TEST(CoordinateUtils, NormalizeAngle)
{
  EXPECT_NEAR(normalize_angle(0.0), 0.0, 1e-10);
  EXPECT_NEAR(normalize_angle(M_PI), M_PI, 1e-10);
  EXPECT_NEAR(normalize_angle(-M_PI), M_PI, 1e-10);
  EXPECT_NEAR(normalize_angle(2.0 * M_PI), 0.0, 1e-10);
  EXPECT_NEAR(normalize_angle(-2.0 * M_PI), 0.0, 1e-10);
  EXPECT_NEAR(normalize_angle(M_PI / 2.0), M_PI / 2.0, 1e-10);
  EXPECT_NEAR(normalize_angle(-M_PI / 2.0), -M_PI / 2.0, 1e-10);
}

TEST(CoordinateUtils, NormalizeAnglePositive)
{
  EXPECT_NEAR(normalize_angle_positive(0.0), 0.0, 1e-10);
  EXPECT_NEAR(normalize_angle_positive(M_PI), M_PI, 1e-10);
  EXPECT_NEAR(normalize_angle_positive(-M_PI), M_PI, 1e-10);
  EXPECT_NEAR(normalize_angle_positive(2.0 * M_PI), 0.0, 1e-10);
  EXPECT_NEAR(normalize_angle_positive(M_PI / 2.0), M_PI / 2.0, 1e-10);
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
