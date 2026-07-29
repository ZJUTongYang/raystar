#include <gtest/gtest.h>

#include <raystar/polymap.h>

#include <algorithm>
#include <iterator>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace raystar;

namespace
{

using FT = exact_geometry::FT;
using Point = exact_geometry::Point;

Polymap makeReadyPolymap(
  const GridMap& map, int start_x, int start_y, int goal_x, int goal_y)
{
  auto created = Polymap::create(map, start_x, start_y, goal_x, goal_y);
  if (!created)
  {
    throw std::runtime_error(
      "Polymap construction failed with status " +
      std::to_string(static_cast<int>(created.status)) + ": " + created.error);
  }
  return std::move(*created.value);
}

struct Direction
{
  FT x;
  FT y;
};

struct Segment
{
  Point from;
  Point to;
};

struct RaySegmentHit
{
  FT ray_parameter;
  FT segment_parameter;
};

struct RayCastResult
{
  bool collinear = false;
  bool endpoint_hit = false;
  std::vector<FT> positive_parameters;
};

enum class DomainLocation
{
  blocked,
  free,
  boundary
};

FT cross(const Direction& first, const Direction& second)
{
  return first.x * second.y - first.y * second.x;
}

FT dot(const Direction& first, const Direction& second)
{
  return first.x * second.x + first.y * second.y;
}

FT exactAbs(const FT& value)
{
  return value < FT(0) ? -value : value;
}

Direction directionFrom(const Point& source, const Point& target)
{
  return {target.x() - source.x(), target.y() - source.y()};
}

FT squaredDistance(const Point& first, const Point& second)
{
  const Direction difference = directionFrom(first, second);
  return dot(difference, difference);
}

bool isZero(const Direction& direction)
{
  return direction.x == FT(0) && direction.y == FT(0);
}

int angularHalf(const Direction& direction)
{
  return direction.y > FT(0) ||
         (direction.y == FT(0) && direction.x >= FT(0)) ? 0 : 1;
}

bool angleLess(const Direction& first, const Direction& second)
{
  const int first_half = angularHalf(first);
  const int second_half = angularHalf(second);
  if (first_half != second_half)
    return first_half < second_half;

  const FT determinant = cross(first, second);
  if (determinant != FT(0))
    return determinant > FT(0);
  return dot(first, first) < dot(second, second);
}

bool sameRay(const Direction& first, const Direction& second)
{
  return cross(first, second) == FT(0) && dot(first, second) > FT(0);
}

Direction l1Normalized(const Direction& direction)
{
  const FT length = exactAbs(direction.x) + exactAbs(direction.y);
  return {direction.x / length, direction.y / length};
}

std::optional<RaySegmentHit> intersectRaySegment(
  const Point& source, const Direction& direction, const Segment& segment,
  bool& collinear)
{
  collinear = false;
  const Direction edge{
    segment.to.x() - segment.from.x(),
    segment.to.y() - segment.from.y()};
  if (isZero(edge))
    return std::nullopt;

  const Direction offset{
    segment.from.x() - source.x(),
    segment.from.y() - source.y()};
  const FT determinant = cross(direction, edge);
  if (determinant == FT(0))
  {
    if (cross(offset, direction) == FT(0))
    {
      const FT squared_length = dot(direction, direction);
      const Direction to_offset{
        segment.to.x() - source.x(),
        segment.to.y() - source.y()};
      const FT from_parameter = dot(offset, direction) / squared_length;
      const FT to_parameter = dot(to_offset, direction) / squared_length;
      collinear = from_parameter > FT(0) || to_parameter > FT(0);
    }
    return std::nullopt;
  }

  const FT ray_parameter = cross(offset, edge) / determinant;
  const FT segment_parameter = cross(offset, direction) / determinant;
  if (ray_parameter < FT(0) ||
      segment_parameter < FT(0) || segment_parameter > FT(1))
  {
    return std::nullopt;
  }
  return RaySegmentHit{ray_parameter, segment_parameter};
}

RayCastResult castRay(
  const Point& source, const Direction& direction,
  const std::vector<Segment>& segments)
{
  RayCastResult result;
  for (const auto& segment : segments)
  {
    bool collinear = false;
    const auto hit = intersectRaySegment(
      source, direction, segment, collinear);
    result.collinear = result.collinear || collinear;
    if (!hit || hit->ray_parameter <= FT(0))
      continue;
    if (hit->segment_parameter == FT(0) ||
        hit->segment_parameter == FT(1))
    {
      result.endpoint_hit = true;
    }
    result.positive_parameters.emplace_back(hit->ray_parameter);
  }

  std::sort(result.positive_parameters.begin(),
    result.positive_parameters.end());
  result.positive_parameters.erase(std::unique(
    result.positive_parameters.begin(), result.positive_parameters.end()),
    result.positive_parameters.end());
  return result;
}

std::vector<Segment> obstacleSegments(const Polymap& polymap)
{
  std::vector<Segment> segments;
  for (const auto& obstacle : polymap.obstacles())
  {
    const auto& vertices = obstacle.ordered_vertices_;
    for (size_t index = 0; index < vertices.size(); ++index)
    {
      const auto& from = vertices[index];
      const auto& to = vertices[(index + 1) % vertices.size()];
      segments.push_back(Segment{
        Point(from.first, from.second), Point(to.first, to.second)});
    }
  }
  return segments;
}

std::vector<Segment> visibilityBoundarySegments(
  const std::vector<Segment>& obstacle_segments, const Point& source,
  const VisibilityRegion& visibility_region)
{
  std::vector<Segment> segments;
  if (visibility_region.empty())
    return segments;

  const bool open_sector = std::any_of(
    obstacle_segments.begin(), obstacle_segments.end(),
    [&source](const Segment& segment) {
      return segment.from == source || segment.to == source;
    });
  if (open_sector)
    segments.push_back(Segment{source, exactPoint(visibility_region.front())});

  for (size_t index = 1; index < visibility_region.size(); ++index)
  {
    segments.push_back(Segment{
      exactPoint(visibility_region[index - 1]),
      exactPoint(visibility_region[index])});
  }

  if (open_sector)
  {
    segments.push_back(Segment{exactPoint(visibility_region.back()), source});
  }
  else
  {
    segments.push_back(Segment{
      exactPoint(visibility_region.back()),
      exactPoint(visibility_region.front())});
  }
  return segments;
}

DomainLocation locateInFreeDomain(
  const Point& query, const std::vector<Segment>& obstacle_segments)
{
  bool inside = false;
  for (const auto& segment : obstacle_segments)
  {
    const auto exact_segment =
      exact_geometry::Kernel::Segment_2(segment.from, segment.to);
    if (exact_segment.has_on(query))
      return DomainLocation::boundary;

    const bool from_above = segment.from.y() > query.y();
    const bool to_above = segment.to.y() > query.y();
    if (from_above == to_above)
      continue;

    const FT intersection_x = segment.from.x() +
      (segment.to.x() - segment.from.x()) *
      (query.y() - segment.from.y()) /
      (segment.to.y() - segment.from.y());
    if (intersection_x == query.x())
      return DomainLocation::boundary;
    if (intersection_x > query.x())
      inside = !inside;
  }
  return inside ? DomainLocation::free : DomainLocation::blocked;
}

std::vector<Direction> eventDirections(
  const Point& source, const std::vector<Segment>& obstacle_segments,
  const VisibilityRegion& visibility_region)
{
  std::vector<Direction> directions = {
    {FT(1), FT(0)}, {FT(0), FT(1)},
    {FT(-1), FT(0)}, {FT(0), FT(-1)}};
  for (const auto& segment : obstacle_segments)
  {
    directions.emplace_back(directionFrom(source, segment.from));
    directions.emplace_back(directionFrom(source, segment.to));
  }
  for (const auto& endpoint : visibility_region)
    directions.emplace_back(directionFrom(source, exactPoint(endpoint)));

  directions.erase(std::remove_if(directions.begin(), directions.end(),
    [](const Direction& direction) { return isZero(direction); }),
    directions.end());
  std::sort(directions.begin(), directions.end(), angleLess);

  std::vector<Direction> unique_directions;
  for (const auto& direction : directions)
  {
    if (unique_directions.empty() ||
        !sameRay(unique_directions.back(), direction))
    {
      unique_directions.emplace_back(direction);
    }
  }
  return unique_directions;
}

std::vector<Direction> sectorProbeDirections(
  const std::vector<Direction>& events)
{
  std::vector<Direction> probes;
  probes.reserve(events.size() * 2);
  for (size_t index = 0; index < events.size(); ++index)
  {
    const Direction first = l1Normalized(events[index]);
    const Direction second =
      l1Normalized(events[(index + 1) % events.size()]);
    probes.push_back(Direction{
      FT(2) * first.x + second.x,
      FT(2) * first.y + second.y});
    probes.push_back(Direction{
      first.x + FT(2) * second.x,
      first.y + FT(2) * second.y});
  }
  return probes;
}

void expectVisibilityMatchesIndependentRayOracle(
  Polymap& polymap, int source_x, int source_y)
{
  VisibilityRegion visibility_region;
  std::string error;
  ASSERT_TRUE(polymap.getVisibilityRegion(
    source_x, source_y, visibility_region, &error)) << error;
  ASSERT_GE(visibility_region.size(), 2u);

  const Point source(source_x, source_y);
  const auto obstacle_segments = obstacleSegments(polymap);
  const auto boundary_segments = visibilityBoundarySegments(
    obstacle_segments, source, visibility_region);
  const auto events = eventDirections(
    source, obstacle_segments, visibility_region);
  ASSERT_GE(events.size(), 4u);
  const auto probes = sectorProbeDirections(events);
  ASSERT_GE(probes.size(), 8u);

  size_t checked_probe_count = 0;
  for (const auto& direction : probes)
  {
    SCOPED_TRACE(testing::Message()
      << "source=(" << source_x << "," << source_y << ")"
      << " direction=(" << CGAL::to_double(direction.x) << ","
      << CGAL::to_double(direction.y) << ")");

    const auto obstacle_cast = castRay(
      source, direction, obstacle_segments);
    ASSERT_FALSE(obstacle_cast.collinear)
      << "Probe generation must avoid obstacle-edge event rays";
    ASSERT_FALSE(obstacle_cast.endpoint_hit)
      << "Probe generation must avoid obstacle vertices";
    ASSERT_FALSE(obstacle_cast.positive_parameters.empty())
      << "The bordered map must stop every ray";

    const FT reference_extent = obstacle_cast.positive_parameters.front();
    const Point midpoint(
      source.x() + direction.x * reference_extent / FT(2),
      source.y() + direction.y * reference_extent / FT(2));
    const DomainLocation midpoint_location = locateInFreeDomain(
      midpoint, obstacle_segments);
    ASSERT_NE(midpoint_location, DomainLocation::boundary);
    const bool reference_visible = midpoint_location == DomainLocation::free;

    const auto boundary_cast = castRay(
      source, direction, boundary_segments);
    ASSERT_FALSE(boundary_cast.collinear)
      << "Probe generation must avoid visibility-boundary event rays";
    ASSERT_FALSE(boundary_cast.endpoint_hit)
      << "Probe generation must avoid visibility endpoints";

    if (!reference_visible)
    {
      EXPECT_TRUE(boundary_cast.positive_parameters.empty())
        << "The open visibility fan must not cover the blocked source wedge";
    }
    else
    {
      ASSERT_EQ(boundary_cast.positive_parameters.size(), 1u)
        << "A source-visible ray must cross the star-shaped boundary once";
      EXPECT_TRUE(boundary_cast.positive_parameters.front() == reference_extent)
        << "candidate="
        << CGAL::to_double(boundary_cast.positive_parameters.front())
        << " reference=" << CGAL::to_double(reference_extent);
    }
    ++checked_probe_count;
  }
  EXPECT_EQ(checked_probe_count, probes.size());
}

GridMap makeBorderedMap(unsigned int width, unsigned int height)
{
  GridMap map;
  map.width = width;
  map.height = height;
  map.resolution = 1.0f;
  map.data.assign(static_cast<size_t>(width) * height, 0);
  for (unsigned int x = 0; x < width; ++x)
  {
    map.data[x] = 1;
    map.data[(height - 1) * width + x] = 1;
  }
  for (unsigned int y = 0; y < height; ++y)
  {
    map.data[y * width] = 1;
    map.data[y * width + width - 1] = 1;
  }
  return map;
}

void fillRectangle(
  GridMap& map, unsigned int min_x, unsigned int min_y,
  unsigned int max_x, unsigned int max_y)
{
  for (unsigned int y = min_y; y < max_y; ++y)
    for (unsigned int x = min_x; x < max_x; ++x)
      map.data[y * map.width + x] = 1;
}

}  // namespace

TEST(VisibilityReference, ClosedInteriorSourcesMatchNaiveExactRayCasting)
{
  auto empty_map = makeBorderedMap(12, 12);
  Polymap empty_polymap = makeReadyPolymap(empty_map, 3, 4, 9, 8);
  ASSERT_TRUE(empty_polymap.hasSolution());
  ASSERT_TRUE(empty_polymap.isCDTReady());
  expectVisibilityMatchesIndependentRayOracle(empty_polymap, 3, 4);

  auto rectangle_map = makeBorderedMap(16, 14);
  fillRectangle(rectangle_map, 6, 4, 9, 8);
  Polymap rectangle_polymap = makeReadyPolymap(rectangle_map, 2, 6, 13, 10);
  ASSERT_TRUE(rectangle_polymap.hasSolution());
  ASSERT_TRUE(rectangle_polymap.isCDTReady());
  expectVisibilityMatchesIndependentRayOracle(rectangle_polymap, 2, 6);

  // Obstacle simplification preserves the constructor start/goal, so build a
  // separate instance when exercising a source on the other side of the map.
  Polymap rectangle_from_right = makeReadyPolymap(rectangle_map, 12, 3, 2, 6);
  ASSERT_TRUE(rectangle_from_right.hasSolution());
  ASSERT_TRUE(rectangle_from_right.isCDTReady());
  expectVisibilityMatchesIndependentRayOracle(rectangle_from_right, 12, 3);

  auto concave_map = makeBorderedMap(18, 16);
  fillRectangle(concave_map, 7, 5, 9, 12);
  fillRectangle(concave_map, 7, 5, 13, 7);
  Polymap concave_polymap = makeReadyPolymap(concave_map, 3, 8, 11, 9);
  ASSERT_TRUE(concave_polymap.hasSolution());
  ASSERT_TRUE(concave_polymap.isCDTReady());
  expectVisibilityMatchesIndependentRayOracle(concave_polymap, 3, 8);
}

TEST(VisibilityReference, ObstacleVertexOpenFansMatchNaiveExactRayCasting)
{
  auto rectangle_map = makeBorderedMap(16, 14);
  fillRectangle(rectangle_map, 6, 4, 9, 8);
  Polymap rectangle_polymap = makeReadyPolymap(rectangle_map, 2, 6, 13, 10);
  ASSERT_TRUE(rectangle_polymap.isCDTReady());
  ASSERT_TRUE(rectangle_polymap.isValidTopology(
    rectangle_polymap.locateVertex(6, 4)));
  expectVisibilityMatchesIndependentRayOracle(rectangle_polymap, 6, 4);

  auto concave_map = makeBorderedMap(18, 16);
  fillRectangle(concave_map, 7, 5, 9, 12);
  fillRectangle(concave_map, 7, 5, 13, 7);
  Polymap concave_polymap = makeReadyPolymap(concave_map, 3, 8, 11, 9);
  ASSERT_TRUE(concave_polymap.isCDTReady());
  ASSERT_TRUE(concave_polymap.isValidTopology(
    concave_polymap.locateVertex(9, 10)));
  expectVisibilityMatchesIndependentRayOracle(concave_polymap, 9, 10);
}

TEST(VisibilityReference, ClosedSameRayDiscontinuitiesKeepDirectedOrder)
{
  auto rectangle_map = makeBorderedMap(16, 14);
  fillRectangle(rectangle_map, 6, 4, 9, 8);
  Polymap polymap = makeReadyPolymap(rectangle_map, 2, 6, 13, 10);
  ASSERT_TRUE(polymap.isCDTReady());

  VisibilityRegion visibility_region;
  std::string error;
  ASSERT_TRUE(polymap.getVisibilityRegion(
    2, 6, visibility_region, &error)) << error;

  const Point source(2, 6);
  const auto lower = std::find_if(
    visibility_region.begin(), visibility_region.end(),
    [](const BoundaryEndpoint& endpoint) {
      return exactPoint(endpoint) == Point(6, 4);
    });
  ASSERT_NE(lower, visibility_region.end());
  ASSERT_NE(lower, visibility_region.begin());
  const Point lower_near = exactPoint(*lower);
  const Point lower_far = exactPoint(*std::prev(lower));
  EXPECT_TRUE(exact_geometry::isSameDirectedRay(
    source, lower_far, lower_near));
  EXPECT_TRUE(squaredDistance(source, lower_far) >
    squaredDistance(source, lower_near));

  const auto upper = std::find_if(
    visibility_region.begin(), visibility_region.end(),
    [](const BoundaryEndpoint& endpoint) {
      return exactPoint(endpoint) == Point(6, 8);
    });
  ASSERT_NE(upper, visibility_region.end());
  ASSERT_NE(std::next(upper), visibility_region.end());
  const Point upper_near = exactPoint(*upper);
  const Point upper_far = exactPoint(*std::next(upper));
  EXPECT_TRUE(exact_geometry::isSameDirectedRay(
    source, upper_near, upper_far));
  EXPECT_TRUE(squaredDistance(source, upper_near) <
    squaredDistance(source, upper_far));
}
