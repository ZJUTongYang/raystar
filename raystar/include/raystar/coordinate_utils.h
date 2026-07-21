#pragma once

#include <cmath>
#include <cstdint>
#include <utility>

namespace raystar
{

/// Grid map data structure (binary occupancy: 0=free, 1=occupied)
struct GridMap
{
  std::vector<uint8_t> data;
  unsigned int width = 0;
  unsigned int height = 0;
  float resolution = 0.05f;
  double origin_x = 0.0;
  double origin_y = 0.0;

  /// Access cell value (row-major). Returns 1 (occupied) for out-of-bounds.
  inline uint8_t operator()(unsigned int x, unsigned int y) const
  {
    if (x >= width || y >= height) return 1;
    return data[y * width + x];
  }

  inline uint8_t& operator()(unsigned int x, unsigned int y)
  {
    return data[y * width + x];
  }

  inline uint8_t at(unsigned int x, unsigned int y) const
  {
    if (x >= width || y >= height) return 1;
    return data[y * width + x];
  }
};

/// World coordinates (meters) -> Grid coordinates (cells).
/// Returns false if outside the map.
inline bool worldToMap(const GridMap& map, double wx, double wy,
                       unsigned int& mx, unsigned int& my)
{
  double dx = (wx - map.origin_x) / map.resolution;
  double dy = (wy - map.origin_y) / map.resolution;
  if (dx < 0.0 || dy < 0.0) return false;
  mx = static_cast<unsigned int>(dx);
  my = static_cast<unsigned int>(dy);
  return mx < map.width && my < map.height;
}

/// Grid coordinates (cells) -> World coordinates (meters).
inline void mapToWorld(const GridMap& map, unsigned int mx, unsigned int my,
                       double& wx, double& wy)
{
  wx = map.origin_x + static_cast<double>(mx) * map.resolution;
  wy = map.origin_y + static_cast<double>(my) * map.resolution;
}

/// Normalize angle to [-PI, PI]
inline double normalize_angle(double angle)
{
  while (angle > M_PI) angle -= 2.0 * M_PI;
  while (angle <= -M_PI) angle += 2.0 * M_PI;
  return angle;
}

/// Normalize angle to [0, 2*PI)
inline double normalize_angle_positive(double angle)
{
  while (angle >= 2.0 * M_PI) angle -= 2.0 * M_PI;
  while (angle < 0.0) angle += 2.0 * M_PI;
  return angle;
}

}  // namespace raystar
