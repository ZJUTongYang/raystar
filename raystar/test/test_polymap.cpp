#include <gtest/gtest.h>
#include <raystar/polymap.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <set>
#include <stack>
#include <utility>
#include <vector>
#include <CGAL/exceptions.h>

namespace raystar {

class PolymapTestPeer {
public:
  using CDT = constrained_delaunay_triangulation::CDT;
  using Validator = bool (*)(const CDT&);

  static bool reconstructWithValidator(Polymap& polymap, Validator validator) {
    return polymap.constructCGALRelated(validator, polymap.construction_error_);
  }

  static const std::string& constructionError(const Polymap& polymap) {
    return polymap.construction_error_;
  }

  static bool constructionStopped(const Polymap& polymap) {
    return polymap.construction_stopped_;
  }

  static std::vector<Obs>& mutableObstacles(Polymap& polymap) {
    return polymap.obstacles();
  }

  static bool refreshObstacles(Polymap& polymap, int start_x, int start_y, int goal_x, int goal_y) {
    return polymap.getPolyObstacles(start_x, start_y, goal_x, goal_y);
  }

  static Polymap constructUnchecked(
    const GridMap& map, int start_x, int start_y, int goal_x, int goal_y) {
    return Polymap(map, start_x, start_y, goal_x, goal_y, StopToken{});
  }

  static size_t cdtVertexCount(const Polymap& polymap) {
    return polymap.cdt_.number_of_vertices();
  }

  static size_t cdtTableSize(const Polymap& polymap) {
    return polymap.cdt_table_.size();
  }

  static size_t facetCount(const Polymap& polymap) {
    return polymap.facets_.size();
  }

  static size_t visibilityCacheSize(const Polymap& polymap) {
    return polymap.visibility_storage_.size();
  }

  static std::pair<size_t, size_t> registrySizes(const Polymap& polymap) {
    return {polymap.vertices_location_x_flat_.size(), polymap.vertices_location_y_flat_.size()};
  }

  static void clearVisibilityCache(Polymap& polymap) {
    polymap.visibility_storage_.clear();
  }

  static bool reverseCachedVisibility(Polymap& polymap, int source_x, int source_y) {
    const int key = source_x + source_y * polymap.width();
    const auto cached = polymap.visibility_storage_.find(key);
    if (cached == polymap.visibility_storage_.end())
      return false;
    std::reverse(cached->second.vertices.begin(), cached->second.vertices.end());
    return true;
  }

  static int recordedCdtVertexCount(const Polymap& polymap) {
    return polymap.cdt_ver_num_;
  }

  static bool validateObstacleTopology(const Polymap& polymap, std::string& error) {
    return polymap.validateObstacleTopology(error, true);
  }

  static bool simplificationChordIsTopologicallySafe(const Polymap& polymap,
                                                     size_t obstacle_index,
                                                     size_t current_index) {
    return polymap.simplificationChordIsTopologicallySafe(obstacle_index, current_index);
  }
};

}  // namespace raystar

using namespace raystar;

namespace {

template <typename... Args>
Polymap makeReadyPolymap(Args&&... args) {
  auto result = Polymap::create(std::forward<Args>(args)...);
  if (!result) {
    throw std::runtime_error(result.error.empty() ? "Polymap construction did not produce a value"
                                                  : result.error);
  }
  return std::move(*result.value);
}

double testPathLength(const std::vector<Point2d>& path) {
  double result = 0.0;
  for (size_t index = 1; index < path.size(); ++index) {
    result += std::hypot(path[index].first - path[index - 1].first,
                         path[index].second - path[index - 1].second);
  }
  return result;
}

bool rejectCDT(const PolymapTestPeer::CDT&) {
  return false;
}

bool throwFromCDTValidator(const PolymapTestPeer::CDT&) {
  throw CGAL::Assertion_exception(
    "CGAL", "injected validator assertion", __FILE__, __LINE__, "injected validator failure");
}

void expectClearedCDTState(const Polymap& polymap) {
  EXPECT_FALSE(polymap.isCDTReady());
  EXPECT_EQ(PolymapTestPeer::cdtVertexCount(polymap), 0u);
  EXPECT_EQ(PolymapTestPeer::recordedCdtVertexCount(polymap), 0);
  EXPECT_EQ(PolymapTestPeer::cdtTableSize(polymap), 0u);
  EXPECT_EQ(PolymapTestPeer::facetCount(polymap), 0u);
  EXPECT_EQ(PolymapTestPeer::visibilityCacheSize(polymap), 0u);
  EXPECT_TRUE(polymap.getCDTEdges().empty());
}

using GridVertex = std::pair<int, int>;
using GridBoundaryEdge = std::pair<GridVertex, GridVertex>;

GridBoundaryEdge makeGridBoundaryEdge(const GridVertex& first, const GridVertex& second) {
  if (first <= second)
    return {first, second};
  return {second, first};
}

struct ExtractedContourEdges {
  std::set<GridBoundaryEdge> edges;
  size_t segment_count = 0;
  bool all_segments_are_unit_length = true;
};

ExtractedContourEdges collectContourEdges(const std::vector<Obs>& obstacles) {
  ExtractedContourEdges result;
  for (const auto& obstacle : obstacles) {
    const auto& vertices = obstacle.ordered_vertices_;
    if (vertices.empty())
      continue;
    for (size_t index = 0; index < vertices.size(); ++index) {
      const auto& from = vertices[index];
      const auto& to = vertices[(index + 1) % vertices.size()];
      result.all_segments_are_unit_length =
        result.all_segments_are_unit_length &&
        std::abs(from.first - to.first) + std::abs(from.second - to.second) == 1;
      result.edges.emplace(makeGridBoundaryEdge(from, to));
      ++result.segment_count;
    }
  }
  return result;
}

long long twiceSignedContourArea(const Obs& obstacle) {
  long long area = 0;
  const auto& vertices = obstacle.ordered_vertices_;
  for (size_t index = 0; index < vertices.size(); ++index) {
    const auto& current = vertices[index];
    const auto& next = vertices[(index + 1) % vertices.size()];
    area += static_cast<long long>(current.first) * next.second -
            static_cast<long long>(next.first) * current.second;
  }
  return area;
}

std::set<GridBoundaryEdge> collectExpectedReachableBoundaryEdges(const GridMap& map,
                                                                 int start_x,
                                                                 int start_y) {
  const int width = static_cast<int>(map.width);
  const int height = static_cast<int>(map.height);
  std::vector<uint8_t> visited(static_cast<size_t>(width) * static_cast<size_t>(height), 0);
  std::stack<GridVertex> pending;
  pending.push({start_x, start_y});
  std::set<GridBoundaryEdge> result;

  while (!pending.empty()) {
    const auto [x, y] = pending.top();
    pending.pop();
    if (x < 0 || y < 0 || x >= width || y >= height)
      continue;
    const int index = x + y * width;
    if (visited[static_cast<size_t>(index)] != 0 || map.data[static_cast<size_t>(index)] != 0) {
      continue;
    }
    visited[static_cast<size_t>(index)] = 1;

    const auto visit_or_record =
      [&](int neighbor_x, int neighbor_y, const GridVertex& edge_from, const GridVertex& edge_to) {
        const int neighbor_index = neighbor_x + neighbor_y * width;
        if (map.data[static_cast<size_t>(neighbor_index)] != 0)
          result.emplace(makeGridBoundaryEdge(edge_from, edge_to));
        else
          pending.push({neighbor_x, neighbor_y});
      };

    if (x > 0)
      visit_or_record(x - 1, y, {x, y}, {x, y + 1});
    if (x + 1 < width)
      visit_or_record(x + 1, y, {x + 1, y}, {x + 1, y + 1});
    if (y > 0)
      visit_or_record(x, y - 1, {x, y}, {x + 1, y});
    if (y + 1 < height)
      visit_or_record(x, y + 1, {x, y + 1}, {x + 1, y + 1});
  }
  return result;
}

GridMap makeBorderedFreeMap(unsigned int width, unsigned int height) {
  GridMap map;
  map.width = width;
  map.height = height;
  map.resolution = 1.0f;
  map.data.assign(static_cast<size_t>(width) * height, 0);
  for (unsigned int x = 0; x < width; ++x) {
    map.data[x] = 1;
    map.data[(height - 1) * width + x] = 1;
  }
  for (unsigned int y = 0; y < height; ++y) {
    map.data[y * width] = 1;
    map.data[y * width + width - 1] = 1;
  }
  return map;
}

}  // namespace

TEST(VisibilityPredicates, ClassifiesAllBungiuCasesAndDegenerateRays) {
  using exact_geometry::Point;
  using exact_geometry::PortalRayPosition;
  const Point source(0, 0);
  const Point lower(10, 0);
  const Point upper(0, 10);

  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, upper, Point(10, -1)),
            PortalRayPosition::before_lower);
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, upper, Point(5, 0)),
            PortalRayPosition::equal_lower);
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, upper, Point(5, 5)),
            PortalRayPosition::inside);
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, upper, Point(0, 5)),
            PortalRayPosition::equal_upper);
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, upper, Point(-1, 10)),
            PortalRayPosition::after_upper);

  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, upper, Point(-5, 0)),
            PortalRayPosition::after_upper)
    << "The ray opposite lower is +pi, not equal to lower";
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, upper, Point(0, -5)),
            PortalRayPosition::before_lower)
    << "The ray opposite upper is clockwise from lower";
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, upper, source),
            PortalRayPosition::invalid);
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, source, upper, Point(1, 1)),
            PortalRayPosition::invalid);
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, source, Point(1, 1)),
            PortalRayPosition::invalid);

  const Point same_upper(20, 0);
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, same_upper, Point(5, 0)),
            PortalRayPosition::equal_lower);
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, same_upper, Point(1, 1)),
            PortalRayPosition::after_upper);

  const Point opposite_upper(-10, 0);
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, opposite_upper, Point(0, 1)),
            PortalRayPosition::inside);
  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, opposite_upper, opposite_upper),
            PortalRayPosition::equal_upper);

  EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, Point(0, -10), Point(1, 1)),
            PortalRayPosition::invalid)
    << "A portal wider than pi violates the Bungiu CASE representation";
}

TEST(VisibilityPredicates, ExactCaseClassifierMatchesLegacyAwayFromBoundaries) {
  using exact_geometry::Point;
  using exact_geometry::PortalRayPosition;
  const Point source(0, 0);
  std::vector<std::pair<int, int>> directions;
  for (int x = -3; x <= 3; ++x)
    for (int y = -3; y <= 3; ++y)
      if (x != 0 || y != 0)
        directions.emplace_back(x, y);

  const auto legacy_classify = [](const auto& lower, const auto& upper, const auto& candidate) {
    const double lower_angle = std::atan2(lower.second, lower.first);
    const double upper_angle =
      lower_angle + normalize_angle_positive(std::atan2(upper.second, upper.first) - lower_angle);
    const double theta =
      lower_angle + normalize_angle(std::atan2(candidate.second, candidate.first) - lower_angle);
    if (theta < lower_angle)
      return PortalRayPosition::before_lower;
    if (theta == lower_angle)
      return PortalRayPosition::equal_lower;
    if (theta > lower_angle && theta < upper_angle)
      return PortalRayPosition::inside;
    if (theta == upper_angle)
      return PortalRayPosition::equal_upper;
    return PortalRayPosition::after_upper;
  };

  for (const auto& lower_direction : directions) {
    const Point lower(lower_direction.first, lower_direction.second);
    for (const auto& upper_direction : directions) {
      const Point upper(upper_direction.first, upper_direction.second);
      if (!exact_geometry::isPortalSpanAtMostPi(source, lower, upper))
        continue;
      for (const auto& candidate_direction : directions) {
        const Point candidate(candidate_direction.first, candidate_direction.second);
        if (exact_geometry::isSameDirectedRay(source, lower, candidate) ||
            exact_geometry::isSameDirectedRay(source, upper, candidate)) {
          continue;
        }
        SCOPED_TRACE(testing::Message()
                     << "lower=(" << lower_direction.first << "," << lower_direction.second
                     << ") upper=(" << upper_direction.first << "," << upper_direction.second
                     << ") candidate=(" << candidate_direction.first << ","
                     << candidate_direction.second << ")");
        EXPECT_EQ(exact_geometry::classifyPortalRay(source, lower, upper, candidate),
                  legacy_classify(lower_direction, upper_direction, candidate_direction));
      }
    }
  }
}

TEST(VisibilityPredicates, PolarComparatorsPreserveWrapAndSameRayOrder) {
  using exact_geometry::Point;
  const Point source(0, 0);

  std::vector<Point> standard = {
    Point(-1, 0), Point(0, 1), Point(1, 0), Point(0, -1), Point(-10, -1)};
  std::stable_sort(standard.begin(), standard.end(), [&](const Point& first, const Point& second) {
    return exact_geometry::compareRaysLikeAtan2(source, first, second) == CGAL::SMALLER;
  });
  const std::vector<Point> expected_standard = {
    Point(-10, -1), Point(0, -1), Point(1, 0), Point(0, 1), Point(-1, 0)};
  EXPECT_EQ(standard, expected_standard);

  const Point near(2, 2);
  const Point far(10, 10);
  EXPECT_EQ(exact_geometry::compareRaysLikeAtan2(source, near, far), CGAL::EQUAL);
  EXPECT_EQ(exact_geometry::compareRaysLikeAtan2(source, far, near), CGAL::EQUAL);
  std::vector<Point> near_then_far = {near, far};
  std::stable_sort(
    near_then_far.begin(), near_then_far.end(), [&](const Point& first, const Point& second) {
      return exact_geometry::compareRaysLikeAtan2(source, first, second) == CGAL::SMALLER;
    });
  EXPECT_EQ(near_then_far, (std::vector<Point>{near, far}));
  std::vector<Point> far_then_near = {far, near};
  std::stable_sort(
    far_then_near.begin(), far_then_near.end(), [&](const Point& first, const Point& second) {
      return exact_geometry::compareRaysLikeAtan2(source, first, second) == CGAL::SMALLER;
    });
  EXPECT_EQ(far_then_near, (std::vector<Point>{far, near}));

  const Point reference(0, 1);
  std::vector<Point> relative = {Point(1, 0), Point(0, -1), Point(-1, 0), reference};
  std::stable_sort(relative.begin(), relative.end(), [&](const Point& first, const Point& second) {
    return exact_geometry::compareRaysCounterClockwiseFrom(source, reference, first, second) ==
           CGAL::SMALLER;
  });
  EXPECT_EQ(relative, (std::vector<Point>{reference, Point(-1, 0), Point(0, -1), Point(1, 0)}));

  const Point same_reference_ray_near(0, 2);
  const Point same_reference_ray_far(0, 20);
  EXPECT_EQ(exact_geometry::compareRaysCounterClockwiseFrom(
              source, reference, same_reference_ray_near, same_reference_ray_far),
            CGAL::EQUAL);
  EXPECT_EQ(exact_geometry::compareRaysCounterClockwiseFrom(
              source, reference, same_reference_ray_far, same_reference_ray_near),
            CGAL::EQUAL);
}

TEST(VisibilityPredicates, SegmentRayIntersectionUsesExactParameterDomains) {
  using exact_geometry::FT;
  using exact_geometry::Point;

  const auto large_intersection = exact_geometry::intersectSegmentWithRay(
    Point(1000000000, -1000000000), Point(1000000000, 1000000000), Point(0, 0), Point(3, 1));
  ASSERT_TRUE(large_intersection.has_value());
  EXPECT_EQ(large_intersection->x(), FT(1000000000));
  EXPECT_EQ(large_intersection->y() * FT(3), FT(1000000000));

  EXPECT_FALSE(
    exact_geometry::intersectSegmentWithRay(Point(-2, -1), Point(-2, 1), Point(0, 0), Point(1, 0))
      .has_value())
    << "A line hit behind the ray source is not a visibility intersection";
  EXPECT_FALSE(
    exact_geometry::intersectSegmentWithRay(Point(2, 1), Point(2, 2), Point(0, 0), Point(1, 0))
      .has_value())
    << "The infinite blocker line must not extend beyond its segment";
  EXPECT_FALSE(
    exact_geometry::intersectSegmentWithRay(Point(0, 1), Point(3, 1), Point(0, 0), Point(1, 0))
      .has_value());
  EXPECT_FALSE(
    exact_geometry::intersectSegmentWithRay(Point(2, 0), Point(3, 0), Point(0, 0), Point(1, 0))
      .has_value())
    << "A collinear overlap has no unique intersection endpoint";
  EXPECT_FALSE(
    exact_geometry::intersectSegmentWithRay(Point(2, 0), Point(2, 0), Point(0, 0), Point(1, 0))
      .has_value());
  EXPECT_FALSE(
    exact_geometry::intersectSegmentWithRay(Point(2, -1), Point(2, 1), Point(0, 0), Point(0, 0))
      .has_value());
}

TEST(VisibilityPredicates, ClosedCounterClockwiseSweepSupportsWideWrap) {
  using exact_geometry::Point;
  const Point source(0, 0);
  const Point start(-1, 1);
  const Point end(1, 1);

  EXPECT_TRUE(exact_geometry::isClosedCounterClockwiseSweepMember(source, start, end, start));
  EXPECT_TRUE(
    exact_geometry::isClosedCounterClockwiseSweepMember(source, start, end, Point(-1, 0)));
  EXPECT_TRUE(
    exact_geometry::isClosedCounterClockwiseSweepMember(source, start, end, Point(0, -1)));
  EXPECT_TRUE(exact_geometry::isClosedCounterClockwiseSweepMember(source, start, end, Point(1, 0)));
  EXPECT_TRUE(exact_geometry::isClosedCounterClockwiseSweepMember(source, start, end, end));
  EXPECT_FALSE(
    exact_geometry::isClosedCounterClockwiseSweepMember(source, start, end, Point(0, 1)));

  EXPECT_TRUE(
    exact_geometry::isClosedCounterClockwiseSweepMember(source, start, Point(-2, 2), Point(-3, 3)));
  EXPECT_FALSE(
    exact_geometry::isClosedCounterClockwiseSweepMember(source, start, Point(-2, 2), Point(1, 0)));
}

TEST(VisibilityPredicates, CollinearMiddleIsRemovableButSpikeIsNot) {
  using exact_geometry::Point;
  const Point previous(1000000000, 1000000000);
  const Point middle(1500000000, 1500000000);
  const Point next(2000000000, 2000000000);
  EXPECT_TRUE(exact_geometry::isRemovableCollinearMiddle(previous, middle, next));
  EXPECT_TRUE(exact_geometry::isRemovableCollinearMiddle(next, middle, previous));

  const Point spike(2000000000, 2000000000);
  const Point backtrack(1500000000, 1500000000);
  EXPECT_FALSE(exact_geometry::isRemovableCollinearMiddle(previous, spike, backtrack));
  EXPECT_FALSE(
    exact_geometry::isRemovableCollinearMiddle(previous, middle, Point(2000000000, 1999999999)));
}

static GridMap makeMapWithObstacle() {
  GridMap map;
  map.width = 30;
  map.height = 30;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.resize(900, 0);

  for (unsigned int y = 10; y < 20; ++y)
    for (unsigned int x = 10; x < 20; ++x) map.data[y * 30 + x] = 1;

  for (unsigned int x = 0; x < map.width; ++x) {
    map.data[x] = 1;
    map.data[(map.height - 1) * map.width + x] = 1;
  }
  for (unsigned int y = 0; y < map.height; ++y) {
    map.data[y * map.width] = 1;
    map.data[y * map.width + map.width - 1] = 1;
  }

  return map;
}

static GridMap makeMapWithTwoObstacles() {
  GridMap map;
  map.width = 40;
  map.height = 30;
  map.resolution = 1.0f;
  map.data.assign(static_cast<size_t>(map.width) * map.height, 0);

  for (unsigned int y = 10; y < 20; ++y) {
    for (unsigned int x = 10; x < 15; ++x) map.data[y * map.width + x] = 1;
    for (unsigned int x = 24; x < 29; ++x) map.data[y * map.width + x] = 1;
  }
  for (unsigned int x = 0; x < map.width; ++x) {
    map.data[x] = 1;
    map.data[(map.height - 1) * map.width + x] = 1;
  }
  for (unsigned int y = 0; y < map.height; ++y) {
    map.data[y * map.width] = 1;
    map.data[y * map.width + map.width - 1] = 1;
  }
  return map;
}

void expectCertifiedMatchingPortalWitness(const HomotopyShorteningResult& result) {
  ASSERT_TRUE(result) << result.message;
  std::string error;
  EXPECT_TRUE(validateReducedDirectedPortalWitness(result.corridor, &error)) << error;
  error.clear();
  EXPECT_TRUE(validateReducedDirectedPortalWitness(result.output_corridor, &error)) << error;
  error.clear();
  EXPECT_TRUE(sameReducedDirectedPortalWitness(result.corridor, result.output_corridor, &error))
    << error;
  EXPECT_EQ(result.corridor.triangle_occurrences, result.output_corridor.triangle_occurrences);
  EXPECT_EQ(result.corridor.portals.size(), result.output_corridor.portals.size());
}

// Self-contained crop of the obstacle component containing the paper map's
// (101,95)->(85,75) canonical turning edge.  Keeping the original grid
// coordinates is intentional: the regression is about binary floating-point
// interpolation of a non-binary slope at these magnitudes.  The row spans are
// extracted from paper_recovered_from_pdf_9obstacle_cspace_r4.pgm; no runtime
// dependency on the supplementary benchmark package is required.
static GridMap makePaperSlopeBoundaryRegressionMap() {
  GridMap map;
  map.width = 150;
  map.height = 110;
  map.resolution = 1.0f;
  map.data.assign(static_cast<size_t>(map.width) * map.height, 0);
  for (unsigned int x = 0; x < map.width; ++x) {
    map.data[x] = 1;
    map.data[(map.height - 1) * map.width + x] = 1;
  }
  for (unsigned int y = 0; y < map.height; ++y) {
    map.data[y * map.width] = 1;
    map.data[y * map.width + map.width - 1] = 1;
  }

  const std::vector<std::pair<int, int>> occupied_spans{
    {106, 107}, {104, 109}, {103, 110}, {102, 111}, {101, 112}, {101, 112}, {100, 113}, {100, 114},
    {100, 132}, {100, 134}, {100, 135}, {100, 136}, {100, 136}, {100, 137}, {100, 137}, {100, 137},
    {100, 136}, {100, 136}, {97, 135},  {93, 134},  {90, 134},  {88, 133},  {87, 132},  {86, 132},
    {85, 131},  {85, 130},  {84, 129},  {84, 129},  {85, 130},  {85, 131},  {86, 131},  {87, 132},
    {88, 133},  {90, 134},  {93, 134},  {96, 135},  {100, 136}, {100, 136}, {100, 137}, {100, 137},
    {100, 136}, {100, 136}, {100, 135}, {100, 134}, {100, 133}, {100, 131}, {100, 114}, {100, 113},
    {101, 112}, {101, 112}, {102, 111}, {103, 110}, {104, 109}, {106, 107}};
  for (size_t row = 0; row < occupied_spans.size(); ++row) {
    const unsigned int y = static_cast<unsigned int>(45 + row);
    for (int x = occupied_spans[row].first; x <= occupied_spans[row].second; ++x)
      map.data[y * map.width + static_cast<unsigned int>(x)] = 1;
  }
  return map;
}

TEST(Polymap, ConstructionDetectsObstacles) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);

  EXPECT_TRUE(poly.hasSolution());
  EXPECT_TRUE(poly.isCDTReady()) << PolymapTestPeer::constructionError(poly);
  EXPECT_TRUE(PolymapTestPeer::constructionError(poly).empty());
  EXPECT_GE(poly.obstacles().size(), 1u);
  EXPECT_FALSE(poly.getCDTEdges().empty());
}

TEST(Polymap, MultiGoalConstructionProtectsAllReachableEndpoints) {
  const auto map = makeMapWithObstacle();
  const std::vector<PolymapEndpoint> goals = {
    {25, 5, {25.25, 5.25}}, {25, 25, {25.75, 25.75}}, {5, 25, {5.5, 25.5}}};

  auto result = Polymap::create(map, 5, 5, Point2d{5.25, 5.25}, goals, StopToken{});

  ASSERT_EQ(result.status, PolymapCreateStatus::ready) << result.error;
  ASSERT_TRUE(result.value.has_value());
  EXPECT_TRUE(result.value->hasSolution());
  EXPECT_TRUE(result.value->isCDTReady());
  std::string error;
  EXPECT_TRUE(result.value->validateFreeSpaceInterior(Point2d{5.25, 5.25}, &error)) << error;
  for (const auto& goal : goals) {
    error.clear();
    EXPECT_TRUE(result.value->validateFreeSpaceInterior(goal.position, &error)) << error;
  }
}

TEST(Polymap, HomotopyFunnelCrossesTriangleInteriorsInsteadOfFollowingMeshEdges) {
  const auto map = makeBorderedFreeMap(12, 8);
  const Point2d start{2.5, 2.5};
  const Point2d goal{9.5, 5.5};
  auto polymap = makeReadyPolymap(map, 2, 2, 9, 5, start, goal);

  const std::vector<Point2d> reference{start, {2.5, 5.5}, goal};
  const auto result = polymap.shortenPathWithinHomotopy(reference);

  ASSERT_TRUE(result) << result.message;
  ASSERT_EQ(result.path.size(), 2u);
  EXPECT_EQ(result.path.front(), start);
  EXPECT_EQ(result.path.back(), goal);
  EXPECT_NEAR(result.path_cost, std::hypot(7.0, 3.0), 1.0e-12);
  for (const auto& edge : polymap.getCDTEdges()) {
    const Point2d edge_start{static_cast<double>(edge.a.first), static_cast<double>(edge.a.second)};
    const Point2d edge_goal{static_cast<double>(edge.b.first), static_cast<double>(edge.b.second)};
    EXPECT_FALSE((edge_start == start && edge_goal == goal) ||
                 (edge_start == goal && edge_goal == start))
      << "The straight output must not be a triangulation-edge graph path";
  }
  EXPECT_TRUE(result.collision_free);
  EXPECT_TRUE(result.homotopy_preserved);
  EXPECT_TRUE(result.locally_shortest);
  EXPECT_GT(polymap.freeTriangleCount(), 0u);
}

TEST(Polymap, CanonicalPaperSlopeSurvivesStrictTraceWithoutDenseInterpolationDrift) {
  const auto map = makePaperSlopeBoundaryRegressionMap();
  auto polymap = makeReadyPolymap(map, 20, 20, 20, 105, Point2d{20.5, 20.5}, Point2d{20.5, 105.5});
  const Point2d canonical_start{101.0, 95.0};
  const Point2d canonical_goal{85.0, 75.0};
  const std::vector<Point2d> topology_path{canonical_start, canonical_goal};

  const auto obstacle_with_canonical_edge =
    std::find_if(polymap.obstacles().begin(), polymap.obstacles().end(), [](const auto& obstacle) {
      const auto& vertices = obstacle.ordered_vertices_;
      for (size_t index = 0; index < vertices.size(); ++index) {
        if (vertices[index] == std::make_pair(101, 95) &&
            vertices[(index + 1) % vertices.size()] == std::make_pair(85, 75)) {
          return true;
        }
      }
      return false;
    });
  ASSERT_NE(obstacle_with_canonical_edge, polymap.obstacles().end())
    << "The fixture must retain the paper map's exact canonical turning edge";

  const auto canonical = polymap.shortenPathWithinHomotopy(topology_path);
  ASSERT_TRUE(canonical) << canonical.message;
  EXPECT_TRUE(canonical.collision_free);
  EXPECT_TRUE(canonical.homotopy_preserved);
  EXPECT_TRUE(canonical.locally_shortest);

  // Replicate the dense ROS serialization contract: ceil(metric length)
  // equal parameter intervals, followed by the exact endpoint.  The first
  // fractional sample is not exactly collinear in binary64 and falls on the
  // occupied side of this boundary edge under strict CDT tracing.
  const size_t sample_count = static_cast<size_t>(std::ceil(std::hypot(
    canonical_goal.first - canonical_start.first, canonical_goal.second - canonical_start.second)));
  std::vector<Point2d> dense_path;
  dense_path.reserve(sample_count + 1);
  for (size_t sample = 0; sample < sample_count; ++sample) {
    const double fraction = static_cast<double>(sample) / static_cast<double>(sample_count);
    dense_path.emplace_back(
      canonical_start.first + (canonical_goal.first - canonical_start.first) * fraction,
      canonical_start.second + (canonical_goal.second - canonical_start.second) * fraction);
  }
  dense_path.push_back(canonical_goal);
  ASSERT_EQ(sample_count, 26u);
  ASSERT_GE(dense_path.size(), 2u);
  EXPECT_DOUBLE_EQ(dense_path[1].first, 100.38461538461539);
  EXPECT_DOUBLE_EQ(dense_path[1].second, 94.230769230769226);
  const double signed_cross = (canonical_goal.first - canonical_start.first) *
                                (dense_path[1].second - canonical_start.second) -
                              (canonical_goal.second - canonical_start.second) *
                                (dense_path[1].first - canonical_start.first);
  EXPECT_NE(signed_cross, 0.0);

  const auto dense = polymap.shortenPathWithinHomotopy(dense_path);
  EXPECT_EQ(dense.status, HomotopyShorteningStatus::no_corridor);
  EXPECT_NE(dense.message.find("Reference segment 0->1"), std::string::npos);
}

TEST(Polymap, HomotopyFunnelTightensAReferenceBelowAnObstacle) {
  const auto map = makeMapWithObstacle();
  const Point2d start{5.5, 15.5};
  const Point2d goal{25.5, 15.5};
  auto polymap = makeReadyPolymap(map, 5, 15, 25, 15, start, goal);
  const std::vector<Point2d> reference{start, {5.5, 5.5}, {10.0, 10.0}, {20.0, 10.0}, goal};

  const auto result = polymap.shortenPathWithinHomotopy(reference);

  ASSERT_TRUE(result) << result.message;
  ASSERT_GE(result.path.size(), 3u);
  EXPECT_EQ(result.path.front(), start);
  EXPECT_EQ(result.path.back(), goal);
  EXPECT_LT(result.path_cost, testPathLength(reference));
  EXPECT_NE(std::find(result.path.begin(), result.path.end(), Point2d{10.0, 10.0}),
            result.path.end());
  EXPECT_NE(std::find(result.path.begin(), result.path.end(), Point2d{20.0, 10.0}),
            result.path.end());
}

TEST(Polymap, HomotopyFunnelIsCostSymmetricWhenReferenceIsReversed) {
  const auto map = makeMapWithObstacle();
  const Point2d start{5.5, 15.5};
  const Point2d goal{25.5, 15.5};
  auto polymap = makeReadyPolymap(map, 5, 15, 25, 15, start, goal);
  std::vector<Point2d> reference{start, {5.5, 5.5}, {10.0, 10.0}, {20.0, 10.0}, goal};

  const auto forward = polymap.shortenPathWithinHomotopy(reference);
  std::reverse(reference.begin(), reference.end());
  const auto reverse = polymap.shortenPathWithinHomotopy(reference);

  ASSERT_TRUE(forward) << forward.message;
  ASSERT_TRUE(reverse) << reverse.message;
  EXPECT_NEAR(forward.path_cost, reverse.path_cost, 1.0e-12);
  auto reversed_output = reverse.path;
  std::reverse(reversed_output.begin(), reversed_output.end());
  EXPECT_EQ(forward.path, reversed_output);

  auto reversed_faces = reverse.corridor.triangle_occurrences;
  std::reverse(reversed_faces.begin(), reversed_faces.end());
  EXPECT_EQ(forward.corridor.triangle_occurrences, reversed_faces);
  ASSERT_EQ(forward.corridor.portals.size(), reverse.corridor.portals.size());
  for (size_t index = 0; index < forward.corridor.portals.size(); ++index) {
    const auto& first = forward.corridor.portals[index];
    const auto& second = reverse.corridor.portals[reverse.corridor.portals.size() - 1 - index];
    EXPECT_EQ(first.from_triangle, second.to_triangle);
    EXPECT_EQ(first.to_triangle, second.from_triangle);
    EXPECT_EQ(first.left, second.right);
    EXPECT_EQ(first.right, second.left);
  }
}

TEST(Polymap, HomotopyFunnelAcceptsAnIdenticalConfigurationAsZeroTransition) {
  const auto map = makeMapWithObstacle();
  const Point2d endpoint{5.5, 15.5};
  auto polymap = makeReadyPolymap(map, 5, 15, 25, 15, endpoint, Point2d{25.5, 15.5});

  const auto result = polymap.shortenPathWithinHomotopy({endpoint, endpoint});

  expectCertifiedMatchingPortalWitness(result);
  ASSERT_EQ(result.path, std::vector<Point2d>{endpoint});
  EXPECT_DOUBLE_EQ(result.path_cost, 0.0);
  ASSERT_EQ(result.corridor.triangle_occurrences.size(), 1u);
  EXPECT_EQ(result.output_corridor.triangle_occurrences, result.corridor.triangle_occurrences);
  EXPECT_TRUE(result.corridor.portals.empty());
  EXPECT_TRUE(result.collision_free);
  EXPECT_TRUE(result.homotopy_preserved);
  EXPECT_TRUE(result.locally_shortest);
}

TEST(Polymap, HomotopyWitnessPreservesSingleObstacleWindingWithSameSideEndpoints) {
  const auto map = makeMapWithObstacle();
  const Point2d start{5.5, 15.5};
  const Point2d goal{6.5, 15.5};
  auto polymap = makeReadyPolymap(map, 5, 15, 6, 15, start, goal);
  const std::vector<Point2d> clockwise_reference{
    start, {10.0, 10.0}, {20.0, 10.0}, {20.0, 20.0}, {10.0, 20.0}, goal};

  const auto result = polymap.shortenPathWithinHomotopy(clockwise_reference);

  expectCertifiedMatchingPortalWitness(result);
  ASSERT_GT(result.corridor.portals.size(), 1u);
  EXPECT_EQ(result.corridor.triangle_occurrences.front(),
            result.corridor.triangle_occurrences.back());
  EXPECT_GE(std::count(result.corridor.triangle_occurrences.begin(),
                       result.corridor.triangle_occurrences.end(),
                       result.corridor.triangle_occurrences.front()),
            2);
  EXPECT_GT(result.path_cost, std::hypot(goal.first - start.first, goal.second - start.second));
  EXPECT_TRUE(result.collision_free);
  EXPECT_TRUE(result.homotopy_preserved);
  EXPECT_TRUE(result.locally_shortest);

  const auto idempotent = polymap.shortenPathWithinHomotopy(result.path);
  expectCertifiedMatchingPortalWitness(idempotent);
  EXPECT_EQ(idempotent.path, result.path);
  EXPECT_EQ(idempotent.corridor.triangle_occurrences, result.corridor.triangle_occurrences);
}

TEST(Polymap, HomotopyWitnessRetainsRepeatedPortalCycleOccurrences) {
  const auto map = makeMapWithObstacle();
  const Point2d start{5.5, 15.5};
  const Point2d goal{6.5, 15.5};
  auto polymap = makeReadyPolymap(map, 5, 15, 6, 15, start, goal);
  const std::vector<Point2d> one_loop{
    start, {10.0, 10.0}, {20.0, 10.0}, {20.0, 20.0}, {10.0, 20.0}, goal};
  const std::vector<Point2d> two_loops{start,
                                       {10.0, 10.0},
                                       {20.0, 10.0},
                                       {20.0, 20.0},
                                       {10.0, 20.0},
                                       start,
                                       {10.0, 10.0},
                                       {20.0, 10.0},
                                       {20.0, 20.0},
                                       {10.0, 20.0},
                                       goal};

  const auto once = polymap.shortenPathWithinHomotopy(one_loop);
  const auto twice = polymap.shortenPathWithinHomotopy(two_loops);

  expectCertifiedMatchingPortalWitness(once);
  expectCertifiedMatchingPortalWitness(twice);
  ASSERT_FALSE(once.corridor.portals.empty());
  EXPECT_EQ(twice.corridor.portals.size(), 2u * once.corridor.portals.size());
  EXPECT_EQ(twice.corridor.triangle_occurrences.front(),
            twice.corridor.triangle_occurrences.back());
  EXPECT_GE(std::count(twice.corridor.triangle_occurrences.begin(),
                       twice.corridor.triangle_occurrences.end(),
                       twice.corridor.triangle_occurrences.front()),
            3);
  std::string error;
  EXPECT_FALSE(sameReducedDirectedPortalWitness(once.corridor, twice.corridor, &error));
  EXPECT_NE(error.find("counts differ"), std::string::npos) << error;
  EXPECT_GT(twice.path_cost, once.path_cost);
}

TEST(Polymap, HomotopyWitnessPreservesANontrivialTwoObstacleLoop) {
  const auto map = makeMapWithTwoObstacles();
  const Point2d start{5.5, 15.5};
  const Point2d goal{6.5, 15.5};
  auto polymap = makeReadyPolymap(map, 5, 15, 6, 15, start, goal);
  const std::vector<Point2d> reference{start,
                                       {10.0, 10.0},
                                       {15.0, 10.0},
                                       {24.0, 10.0},
                                       {29.0, 10.0},
                                       {29.0, 20.0},
                                       {24.0, 20.0},
                                       {15.0, 20.0},
                                       {10.0, 20.0},
                                       goal};

  const auto result = polymap.shortenPathWithinHomotopy(reference);

  expectCertifiedMatchingPortalWitness(result);
  ASSERT_GT(result.corridor.portals.size(), 2u);
  EXPECT_EQ(result.corridor.triangle_occurrences.front(),
            result.corridor.triangle_occurrences.back());
  EXPECT_LT(result.path.size(), reference.size());
  EXPECT_LE(result.path_cost, testPathLength(reference));
  EXPECT_GT(result.path_cost, std::hypot(goal.first - start.first, goal.second - start.second));
}

TEST(Polymap, ReducedDirectedPortalWitnessRejectsEveryStructuralOrGeometricTamper) {
  const auto map = makeMapWithObstacle();
  const Point2d start{5.5, 15.5};
  const Point2d goal{6.5, 15.5};
  auto polymap = makeReadyPolymap(map, 5, 15, 6, 15, start, goal);
  const auto clockwise = polymap.shortenPathWithinHomotopy(
    {start, {10.0, 10.0}, {20.0, 10.0}, {20.0, 20.0}, {10.0, 20.0}, goal});
  const auto counterclockwise = polymap.shortenPathWithinHomotopy(
    {start, {10.0, 20.0}, {20.0, 20.0}, {20.0, 10.0}, {10.0, 10.0}, goal});
  expectCertifiedMatchingPortalWitness(clockwise);
  expectCertifiedMatchingPortalWitness(counterclockwise);
  ASSERT_FALSE(clockwise.corridor.portals.empty());

  std::string error;
  EXPECT_FALSE(
    sameReducedDirectedPortalWitness(clockwise.corridor, counterclockwise.corridor, &error));
  EXPECT_NE(error.find("differs at ordinal"), std::string::npos) << error;

  auto bad_cardinality = clockwise.output_corridor;
  bad_cardinality.triangle_occurrences.pop_back();
  error.clear();
  EXPECT_FALSE(validateReducedDirectedPortalWitness(bad_cardinality, &error));
  EXPECT_NE(error.find("cardinalities"), std::string::npos) << error;
  error.clear();
  EXPECT_FALSE(sameReducedDirectedPortalWitness(clockwise.corridor, bad_cardinality, &error));
  EXPECT_NE(error.find("Candidate"), std::string::npos) << error;

  auto bad_direction = clockwise.output_corridor;
  std::swap(bad_direction.portals.front().from_triangle, bad_direction.portals.front().to_triangle);
  error.clear();
  EXPECT_FALSE(validateReducedDirectedPortalWitness(bad_direction, &error));
  EXPECT_NE(error.find("does not bind"), std::string::npos) << error;

  auto bad_geometry = clockwise.output_corridor;
  bad_geometry.portals.front().left.first = std::nextafter(bad_geometry.portals.front().left.first,
                                                           std::numeric_limits<double>::infinity());
  error.clear();
  ASSERT_TRUE(validateReducedDirectedPortalWitness(bad_geometry, &error)) << error;
  EXPECT_FALSE(sameReducedDirectedPortalWitness(clockwise.corridor, bad_geometry, &error));
  EXPECT_NE(error.find("geometry differs"), std::string::npos) << error;

  const auto first = clockwise.corridor.portals.front();
  DirectedPortal inverse;
  inverse.from_triangle = first.to_triangle;
  inverse.to_triangle = first.from_triangle;
  inverse.left = first.right;
  inverse.right = first.left;
  TriangleCorridor unreduced;
  unreduced.triangle_occurrences = {first.from_triangle, first.to_triangle, first.from_triangle};
  unreduced.portals = {first, inverse};
  error.clear();
  EXPECT_FALSE(validateReducedDirectedPortalWitness(unreduced, &error));
  EXPECT_NE(error.find("unreduced immediate portal reversal"), std::string::npos) << error;
}

TEST(Polymap, HomotopyFunnelIsIdempotent) {
  const auto map = makeMapWithObstacle();
  const Point2d start{5.5, 15.5};
  const Point2d goal{25.5, 15.5};
  auto polymap = makeReadyPolymap(map, 5, 15, 25, 15, start, goal);
  const std::vector<Point2d> reference{start, {5.5, 5.5}, {10.0, 10.0}, {20.0, 10.0}, goal};

  const auto first = polymap.shortenPathWithinHomotopy(reference);
  ASSERT_TRUE(first) << first.message;
  const auto second = polymap.shortenPathWithinHomotopy(first.path);

  ASSERT_TRUE(second) << second.message;
  EXPECT_EQ(second.path, first.path);
  EXPECT_NEAR(second.path_cost, first.path_cost, 1.0e-12);
}

TEST(Polymap, HomotopyFunnelCancelsAnImmediatePortalBacktrack) {
  const auto map = makeBorderedFreeMap(12, 8);
  const Point2d start{2.5, 2.5};
  const Point2d goal{9.5, 5.5};
  auto polymap = makeReadyPolymap(map, 2, 2, 9, 5, start, goal);

  const auto direct = polymap.shortenPathWithinHomotopy({start, goal});
  const auto with_backtrack = polymap.shortenPathWithinHomotopy({start, {2.5, 5.5}, start, goal});

  ASSERT_TRUE(direct) << direct.message;
  ASSERT_TRUE(with_backtrack) << with_backtrack.message;
  EXPECT_EQ(with_backtrack.path, direct.path);
  EXPECT_EQ(with_backtrack.corridor.triangle_occurrences, direct.corridor.triangle_occurrences);
}

TEST(Polymap, HomotopyFunnelHandlesAReferenceCollinearWithAnInternalPortal) {
  const auto map = makeMapWithObstacle();
  auto polymap = makeReadyPolymap(map, 5, 15, 25, 15, Point2d{5.5, 15.5}, Point2d{25.5, 15.5});
  const auto edge = std::find_if(
    polymap.triangle_edges_.begin(),
    polymap.triangle_edges_.end(),
    [&polymap](const auto& candidate) {
      return !candidate.constrained && candidate.faces[0] >= 0 && candidate.faces[1] >= 0 &&
             polymap.triangle_faces_[static_cast<size_t>(candidate.faces[0])].is_free &&
             polymap.triangle_faces_[static_cast<size_t>(candidate.faces[1])].is_free;
    });
  ASSERT_NE(edge, polymap.triangle_edges_.end());
  const auto interpolate = [&edge](double parameter) {
    return Point2d{edge->a.first + parameter * (edge->b.first - edge->a.first),
                   edge->a.second + parameter * (edge->b.second - edge->a.second)};
  };
  const Point2d start = interpolate(0.25);
  const Point2d middle = interpolate(0.5);
  const Point2d goal = interpolate(0.75);

  const auto result = polymap.shortenPathWithinHomotopy({start, middle, goal});

  ASSERT_TRUE(result) << result.message;
  ASSERT_EQ(result.path, (std::vector<Point2d>{start, goal}));
  EXPECT_NEAR(
    result.path_cost, std::hypot(goal.first - start.first, goal.second - start.second), 1.0e-12);
}

TEST(Polymap, HomotopyFunnelCancellationIsImmediateAndStructured) {
  const auto map = makeBorderedFreeMap(12, 8);
  const Point2d start{2.5, 2.5};
  const Point2d goal{9.5, 5.5};
  auto polymap = makeReadyPolymap(map, 2, 2, 9, 5, start, goal);
  const StopToken stop([]() { return true; });

  const auto result = polymap.shortenPathWithinHomotopy({start, goal}, stop);

  EXPECT_EQ(result.status, HomotopyShorteningStatus::stopped);
  EXPECT_TRUE(result.path.empty());
  EXPECT_FALSE(result.collision_free);
  EXPECT_FALSE(result.homotopy_preserved);
  EXPECT_FALSE(result.locally_shortest);
}

TEST(Polymap, PathTracingReportsStopRequestedInsideSegmentFaceScan) {
  const auto map = makeBorderedFreeMap(12, 8);
  const Point2d start{2.5, 2.5};
  const Point2d goal{9.5, 5.5};
  auto polymap = makeReadyPolymap(map, 2, 2, 9, 5, start, goal);
  size_t poll_count = 0;
  // Poll 1 is the public trace wrapper, poll 2 is the segment loop, and poll
  // 3 is the first triangle-edge iteration inside segment_faces().
  const StopToken stop([&poll_count]() { return ++poll_count == 3; });
  TriangleCorridor corridor;
  corridor.triangle_occurrences.push_back(123u);
  std::string error = "stale diagnostic";

  const auto status = polymap.traceFreeSpacePath({start, goal}, corridor, stop, &error);

  EXPECT_EQ(status, OperationStatus::stopped);
  EXPECT_EQ(poll_count, 3u);
  EXPECT_TRUE(corridor.triangle_occurrences.empty());
  EXPECT_TRUE(corridor.portals.empty());
  EXPECT_TRUE(error.empty());
}

TEST(Polymap, HomotopyFunnelRejectsAReferenceThroughAnObstacle) {
  const auto map = makeMapWithObstacle();
  const Point2d start{5.5, 15.5};
  const Point2d goal{25.5, 15.5};
  auto polymap = makeReadyPolymap(map, 5, 15, 25, 15, start, goal);

  const auto result = polymap.shortenPathWithinHomotopy({start, goal});

  EXPECT_FALSE(result);
  EXPECT_EQ(result.status, HomotopyShorteningStatus::no_corridor);
  EXPECT_FALSE(result.collision_free);
  EXPECT_FALSE(result.homotopy_preserved);
  EXPECT_FALSE(result.locally_shortest);
}

TEST(Polymap, MultiGoalConstructionRejectsAReachabilityMismatchAtomically) {
  auto map = makeBorderedFreeMap(12, 8);
  for (unsigned int y = 1; y + 1 < map.height; ++y) map.data[y * map.width + 6] = 1;
  const std::vector<PolymapEndpoint> goals = {{4, 4, {4.5, 4.5}}, {9, 4, {9.5, 4.5}}};

  auto result = Polymap::create(map, 2, 4, Point2d{2.5, 4.5}, goals, StopToken{});

  EXPECT_EQ(result.status, PolymapCreateStatus::no_path);
  EXPECT_FALSE(result.value.has_value());
  EXPECT_NE(result.error.find("every goal"), std::string::npos);
}

TEST(Polymap, RejectsMalformedGridBeforeCopyingOrSizingRegistries) {
  GridMap malformed;
  malformed.width = 20;
  malformed.height = 20;
  malformed.data.assign(399, 0);

  Polymap integer_endpoints(malformed, 2, 2, 17, 17);
  EXPECT_FALSE(integer_endpoints.hasSolution());
  EXPECT_FALSE(integer_endpoints.isCDTReady());
  EXPECT_FALSE(PolymapTestPeer::constructionStopped(integer_endpoints));
  EXPECT_EQ(PolymapTestPeer::constructionError(integer_endpoints),
            "Invalid map: data size is 399, expected 400");
  EXPECT_EQ(integer_endpoints.width(), 0);
  EXPECT_EQ(integer_endpoints.height(), 0);
  EXPECT_TRUE(integer_endpoints.occupancyData().empty());
  EXPECT_EQ(PolymapTestPeer::registrySizes(integer_endpoints), (std::pair<size_t, size_t>{0, 0}));

  Polymap continuous_endpoints(malformed, 2, 2, 17, 17, Point2d{2.25, 2.25}, Point2d{17.75, 17.75});
  EXPECT_FALSE(continuous_endpoints.hasSolution());
  EXPECT_FALSE(continuous_endpoints.isCDTReady());
  EXPECT_EQ(PolymapTestPeer::constructionError(continuous_endpoints),
            "Invalid map: data size is 399, expected 400");
  EXPECT_TRUE(continuous_endpoints.occupancyData().empty());
  EXPECT_EQ(PolymapTestPeer::registrySizes(continuous_endpoints),
            (std::pair<size_t, size_t>{0, 0}));
}

TEST(Polymap, RejectsUnrepresentableCellCountBeforeRegistryResize) {
  GridMap too_many_cells;
  too_many_cells.width = 46341;
  too_many_cells.height = 46341;

  Polymap polymap(too_many_cells, 0, 0, 0, 0);
  EXPECT_FALSE(polymap.hasSolution());
  EXPECT_FALSE(polymap.isCDTReady());
  EXPECT_EQ(PolymapTestPeer::constructionError(polymap),
            "Invalid map: cell count must fit signed int indexing");
  EXPECT_EQ(polymap.width(), 0);
  EXPECT_EQ(polymap.height(), 0);
  EXPECT_TRUE(polymap.occupancyData().empty());
  EXPECT_EQ(PolymapTestPeer::registrySizes(polymap), (std::pair<size_t, size_t>{0, 0}));
}

TEST(Polymap, CooperativeConstructionStopClearsAllDerivedState) {
  const auto map = makeMapWithObstacle();
  size_t full_poll_count = 0;
  const StopToken count_only([&full_poll_count]() {
    ++full_poll_count;
    return false;
  });
  Polymap complete(map, 5, 5, 25, 25, count_only);
  ASSERT_TRUE(complete.isCDTReady()) << PolymapTestPeer::constructionError(complete);
  ASSERT_GT(full_poll_count, 0u);

  size_t stopped_poll_count = 0;
  const StopToken stop_at_last_poll(
    [&stopped_poll_count, full_poll_count]() { return ++stopped_poll_count >= full_poll_count; });
  Polymap poly(map, 5, 5, 25, 25, stop_at_last_poll);

  EXPECT_EQ(stopped_poll_count, full_poll_count);
  EXPECT_TRUE(PolymapTestPeer::constructionStopped(poly));
  EXPECT_FALSE(poly.hasSolution());
  EXPECT_TRUE(PolymapTestPeer::constructionError(poly).empty());
  EXPECT_TRUE(poly.obstacles().empty());
  EXPECT_EQ(poly.locateVertex(10, 10), std::make_pair(-1, -1));
  expectClearedCDTState(poly);
}

TEST(Polymap, NoSolutionWhenGoalBlocked) {
  GridMap map;
  map.width = 20;
  map.height = 20;
  map.resolution = 1.0f;
  map.origin_x = 0.0;
  map.origin_y = 0.0;
  map.data.resize(400, 0);

  for (unsigned int y = 0; y < 20; ++y) map.data[y * 20 + 10] = 1;

  Polymap poly(map, 5, 10, 15, 10);
  EXPECT_FALSE(poly.hasSolution());
  EXPECT_FALSE(poly.isCDTReady());
  EXPECT_TRUE(PolymapTestPeer::constructionError(poly).empty());
  EXPECT_TRUE(poly.getCDTEdges().empty());
}

TEST(Polymap, FailedObstacleRefreshClearsCommittedDerivedGeometry) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);
  ASSERT_TRUE(poly.hasSolution());
  ASSERT_TRUE(poly.isCDTReady()) << PolymapTestPeer::constructionError(poly);
  ASSERT_NE(poly.locateVertex(10, 10), std::make_pair(-1, -1));

  // The requested goal cell lies inside the central occupied block.  A
  // failed public refresh must not leave a ready CDT/cache referring to an
  // obstacle vector that has just been cleared.
  EXPECT_FALSE(poly.getPolyObstacles(5, 5, 15, 15));
  EXPECT_FALSE(poly.hasSolution());
  EXPECT_FALSE(poly.isCDTReady());
  EXPECT_TRUE(poly.obstacles().empty());
  EXPECT_TRUE(poly.getCDTEdges().empty());
  EXPECT_EQ(poly.locateVertex(10, 10), std::make_pair(-1, -1));

  VisibilityRegion visibility;
  std::string error;
  EXPECT_FALSE(poly.getVisibilityRegion(5, 5, visibility, &error));
  EXPECT_TRUE(visibility.empty());
  EXPECT_FALSE(error.empty());
}

TEST(Polymap, SuccessfulObstacleRefreshInvalidatesOldTopologyAndCDT) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);
  ASSERT_TRUE(poly.hasSolution());
  ASSERT_TRUE(poly.isCDTReady()) << PolymapTestPeer::constructionError(poly);
  ASSERT_NE(poly.locateVertex(10, 10), std::make_pair(-1, -1));

  ASSERT_TRUE(poly.getPolyObstacles(5, 5, 25, 25));
  EXPECT_TRUE(poly.hasSolution());
  EXPECT_FALSE(poly.isCDTReady());
  EXPECT_TRUE(poly.getCDTEdges().empty());
  EXPECT_EQ(poly.locateVertex(10, 10), std::make_pair(-1, -1));

  VisibilityRegion visibility;
  std::string error;
  EXPECT_FALSE(poly.getVisibilityRegion(5, 5, visibility, &error));
  EXPECT_TRUE(visibility.empty());
  EXPECT_NE(error.find("construction failed"), std::string::npos) << error;
}

TEST(Polymap, OddWidthBoundaryEdgesWithEqualLegacySumsRemainDistinct) {
  const auto map = makeBorderedFreeMap(5, 5);
  Polymap polymap(map, 2, 2, 0, 0);
  ASSERT_FALSE(polymap.hasSolution());
  ASSERT_TRUE(polymap.getPolyObstacles(2, 2, 3, 3));
  ASSERT_EQ(polymap.obstacles().size(), 1u);

  const auto extracted = collectContourEdges(polymap.obstacles());
  const auto expected = collectExpectedReachableBoundaryEdges(map, 2, 2);
  EXPECT_TRUE(extracted.all_segments_are_unit_length);
  EXPECT_EQ(extracted.segment_count, 12u);
  EXPECT_EQ(extracted.edges.size(), extracted.segment_count);
  EXPECT_EQ(extracted.edges, expected);
  EXPECT_LT(twiceSignedContourArea(polymap.obstacles().front()), 0)
    << "The outer contour must retain its clockwise free-space-on-right winding";

  // With width 5 these distinct edges had the same old integer key:
  // [6, 11] -> (1,1)-(1,2), [8, 9] -> (3,1)-(4,1), and both summed to 17.
  EXPECT_EQ(extracted.edges.count(makeGridBoundaryEdge({1, 1}, {1, 2})), 1u);
  EXPECT_EQ(extracted.edges.count(makeGridBoundaryEdge({3, 1}, {4, 1})), 1u);
}

TEST(Polymap, ExtractedContoursMatchAllReachableFreeOccupiedInterfaces) {
  auto map = makeBorderedFreeMap(7, 7);
  // This internal contour and the outer contour contain two pairs of edges
  // whose endpoint sums collided under the former width-7 integer key.
  map.data[2 * map.width + 2] = 1;

  Polymap polymap(map, 1, 1, 0, 0);
  ASSERT_FALSE(polymap.hasSolution());
  ASSERT_TRUE(polymap.getPolyObstacles(1, 1, 5, 5));
  ASSERT_EQ(polymap.obstacles().size(), 2u);

  const auto extracted = collectContourEdges(polymap.obstacles());
  const auto expected = collectExpectedReachableBoundaryEdges(map, 1, 1);
  EXPECT_TRUE(extracted.all_segments_are_unit_length);
  EXPECT_EQ(extracted.segment_count, extracted.edges.size())
    << "Each free/occupied interface edge must occur in exactly one contour";
  EXPECT_EQ(extracted.edges, expected);

  size_t clockwise_contours = 0;
  size_t counterclockwise_contours = 0;
  for (const auto& obstacle : polymap.obstacles()) {
    const long long area = twiceSignedContourArea(obstacle);
    clockwise_contours += area < 0 ? 1u : 0u;
    counterclockwise_contours += area > 0 ? 1u : 0u;
  }
  EXPECT_EQ(clockwise_contours, 1u);
  EXPECT_EQ(counterclockwise_contours, 1u);
}

TEST(Polymap, SharedObstacleVertexIsRejectedBeforeRegistryAndCDT) {
  auto map = makeBorderedFreeMap(7, 7);
  map.data[2 * map.width + 2] = 1;
  map.data[3 * map.width + 3] = 1;

  Polymap polymap(map, 4, 2, 5, 5);

  EXPECT_TRUE(polymap.hasSolution())
    << "Grid connectivity must remain independent from topology support";
  EXPECT_FALSE(polymap.isCDTReady());
  EXPECT_TRUE(polymap.getCDTEdges().empty());
  EXPECT_NE(PolymapTestPeer::constructionError(polymap).find("Unsupported obstacle topology:"),
            std::string::npos)
    << PolymapTestPeer::constructionError(polymap);
  EXPECT_NE(PolymapTestPeer::constructionError(polymap).find("(3, 3)"), std::string::npos)
    << PolymapTestPeer::constructionError(polymap);
  EXPECT_NE(PolymapTestPeer::constructionError(polymap).find("obstacle"), std::string::npos)
    << PolymapTestPeer::constructionError(polymap);
  EXPECT_EQ(polymap.locateVertex(3, 3), std::make_pair(-1, -1))
    << "An ambiguous topology must be rejected before registry insertion";
}

TEST(Polymap, RawBoundaryEndpointIsRejectedBeforeTopologySimplification) {
  auto map = makeBorderedFreeMap(7, 7);
  map.data[2 * map.width + 2] = 1;
  map.data[3 * map.width + 3] = 1;

  // The two diagonal obstacle cells create an unsupported shared contour
  // vertex.  Start lies on the raw right edge of cell (2,2), while its floor
  // cell (3,2) is free.  Endpoint validation must win before topology
  // simplification/CGAL so the caller receives the actionable input error.
  Polymap poly(map, 3, 2, 5, 5, Point2d{3.0, 2.5}, Point2d{5.5, 5.5});
  EXPECT_TRUE(poly.hasSolution());
  EXPECT_FALSE(poly.isCDTReady());
  EXPECT_NE(PolymapTestPeer::constructionError(poly).find("Invalid start position"),
            std::string::npos)
    << PolymapTestPeer::constructionError(poly);
  EXPECT_NE(PolymapTestPeer::constructionError(poly).find("boundary"), std::string::npos)
    << PolymapTestPeer::constructionError(poly);
}

TEST(Polymap, TopologyValidatorReportsSharedObstacleIdsAndCoordinate) {
  auto map = makeMapWithObstacle();
  Polymap polymap(map, 5, 5, 25, 25);
  ASSERT_TRUE(polymap.isCDTReady()) << PolymapTestPeer::constructionError(polymap);

  polymap.obstacles() = {Obs{{{0, 0}, {0, 10}, {10, 10}, {10, 0}}},
                         Obs{{{2, 2}, {4, 2}, {4, 4}, {2, 4}}},
                         Obs{{{4, 4}, {6, 4}, {6, 6}, {4, 6}}}};

  std::string error;
  EXPECT_FALSE(PolymapTestPeer::validateObstacleTopology(polymap, error));
  EXPECT_NE(error.find("obstacles 1 and 2 share vertex"), std::string::npos) << error;
  EXPECT_NE(error.find("(4, 4)"), std::string::npos) << error;
}

TEST(Polymap, TopologyValidatorRejectsTjunctionsCrossingsAndOverlaps) {
  auto map = makeMapWithObstacle();
  Polymap polymap(map, 5, 5, 25, 25);
  ASSERT_TRUE(polymap.isCDTReady()) << PolymapTestPeer::constructionError(polymap);

  polymap.obstacles() = {Obs{{{2, 4}, {7, 4}, {7, 8}, {2, 8}}}, Obs{{{3, 2}, {6, 2}, {5, 4}}}};
  std::string error;
  EXPECT_FALSE(PolymapTestPeer::validateObstacleTopology(polymap, error));
  EXPECT_NE(error.find("T-junction"), std::string::npos) << error;

  polymap.obstacles() = {Obs{{{2, 3}, {7, 3}, {7, 8}, {2, 8}}},
                         Obs{{{4, 1}, {9, 1}, {9, 5}, {4, 5}}}};
  error.clear();
  EXPECT_FALSE(PolymapTestPeer::validateObstacleTopology(polymap, error));
  EXPECT_NE(error.find("properly cross"), std::string::npos) << error;

  polymap.obstacles() = {Obs{{{2, 4}, {8, 4}, {8, 8}, {2, 8}}}, Obs{{{4, 4}, {6, 4}, {5, 2}}}};
  error.clear();
  EXPECT_FALSE(PolymapTestPeer::validateObstacleTopology(polymap, error));
  EXPECT_NE(error.find("overlap"), std::string::npos) << error;
}

TEST(Polymap, SimplificationChordRejectsNewIntersection) {
  auto map = makeMapWithObstacle();
  Polymap polymap(map, 5, 5, 25, 25);
  ASSERT_TRUE(polymap.isCDTReady()) << PolymapTestPeer::constructionError(polymap);

  polymap.obstacles() = {Obs{{{2, 2}, {4, 4}, {6, 2}, {6, 6}, {2, 6}}}};
  EXPECT_TRUE(PolymapTestPeer::simplificationChordIsTopologicallySafe(polymap, 0, 1));

  polymap.obstacles().push_back(Obs{{{4, 1}, {5, 3}, {3, 3}}});
  EXPECT_FALSE(PolymapTestPeer::simplificationChordIsTopologicallySafe(polymap, 0, 1));
}

TEST(Polymap, CdtValidatorFalseProducesAtomicConstructionFailure) {
  const auto map = makeMapWithObstacle();
  Polymap polymap(map, 5, 5, 25, 25);

  ASSERT_TRUE(polymap.hasSolution());
  ASSERT_TRUE(polymap.isCDTReady()) << PolymapTestPeer::constructionError(polymap);
  ASSERT_GT(PolymapTestPeer::cdtVertexCount(polymap), 0u);
  ASSERT_GT(PolymapTestPeer::cdtTableSize(polymap), 0u);
  ASSERT_GT(PolymapTestPeer::facetCount(polymap), 0u);
  ASSERT_FALSE(polymap.getCDTEdges().empty());

  VisibilityRegion initial_visibility;
  std::string error;
  ASSERT_TRUE(polymap.getVisibilityRegion(5, 5, initial_visibility, &error)) << error;
  ASSERT_GT(PolymapTestPeer::visibilityCacheSize(polymap), 0u);

  EXPECT_FALSE(PolymapTestPeer::reconstructWithValidator(polymap, &rejectCDT));

  EXPECT_TRUE(polymap.hasSolution())
    << "Grid connectivity and CDT validity must remain independent";
  ASSERT_FALSE(PolymapTestPeer::constructionError(polymap).empty());
  EXPECT_NE(PolymapTestPeer::constructionError(polymap).find("invalid"), std::string::npos);
  expectClearedCDTState(polymap);

  VisibilityRegion visibility = {BoundaryEndpoint{{-1.0, -1.0}, std::monostate{}}};
  error = "stale error";
  EXPECT_FALSE(polymap.getVisibilityRegion(5, 5, visibility, &error));
  EXPECT_TRUE(visibility.empty());
  EXPECT_EQ(error,
            "Visibility region is unavailable because map geometry construction failed: " +
              PolymapTestPeer::constructionError(polymap));

  visibility = {BoundaryEndpoint{{-2.0, -2.0}, std::monostate{}}};
  EXPECT_FALSE(polymap.calculateVisibilityRegion(5, 5, visibility));
  EXPECT_TRUE(visibility.empty());
}

TEST(Polymap, CdtValidatorCGALFailureProducesAtomicConstructionFailure) {
  const auto map = makeMapWithObstacle();
  Polymap polymap(map, 5, 5, 25, 25);

  ASSERT_TRUE(polymap.hasSolution());
  ASSERT_TRUE(polymap.isCDTReady()) << PolymapTestPeer::constructionError(polymap);
  ASSERT_GT(PolymapTestPeer::cdtVertexCount(polymap), 0u);
  ASSERT_FALSE(polymap.getCDTEdges().empty());

  bool reconstruction_succeeded = true;
  EXPECT_NO_THROW(reconstruction_succeeded =
                    PolymapTestPeer::reconstructWithValidator(polymap, &throwFromCDTValidator));

  EXPECT_FALSE(reconstruction_succeeded);
  EXPECT_TRUE(polymap.hasSolution());
  ASSERT_FALSE(PolymapTestPeer::constructionError(polymap).empty());
  EXPECT_NE(PolymapTestPeer::constructionError(polymap).find("injected validator assertion"),
            std::string::npos);
  expectClearedCDTState(polymap);

  VisibilityRegion visibility = {BoundaryEndpoint{{-1.0, -1.0}, std::monostate{}}};
  std::string error;
  EXPECT_FALSE(polymap.getVisibilityRegion(5, 5, visibility, &error));
  EXPECT_TRUE(visibility.empty());
  EXPECT_EQ(error,
            "Visibility region is unavailable because map geometry construction failed: " +
              PolymapTestPeer::constructionError(polymap));
}

TEST(Polymap, VerticesRegistered) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);

  ASSERT_FALSE(poly.obstacles().empty());
  ASSERT_FALSE(poly.obstacles()[0].ordered_vertices_.empty());
  auto first = poly.obstacles()[0].ordered_vertices_[0];
  auto loc = poly.locateVertex(first.first, first.second);
  EXPECT_EQ(loc.first, 0);
  EXPECT_EQ(loc.second, 0);
}

TEST(Polymap, VisibilityRegionComputedAndCachedAtomically) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);

  ASSERT_TRUE(poly.hasSolution());

  VisibilityRegion rich_visibility = {BoundaryEndpoint{{-100.0, -100.0}, std::monostate{}}};
  std::string error = "stale error";
  ASSERT_TRUE(poly.getVisibilityRegion(5, 5, rich_visibility, &error)) << error;
  EXPECT_TRUE(error.empty());
  ASSERT_GE(rich_visibility.size(), 2u);
  ASSERT_TRUE(poly.validateVisibilityRegion(rich_visibility, &error)) << error;
  for (const auto& endpoint : rich_visibility) EXPECT_TRUE(endpoint.exact_position.has_value());

  const auto first_rich_visibility = rich_visibility;

  rich_visibility = {BoundaryEndpoint{{-200.0, -200.0}, std::monostate{}},
                     BoundaryEndpoint{{-201.0, -201.0}, std::monostate{}}};
  error = "another stale error";
  ASSERT_TRUE(poly.getVisibilityRegion(5, 5, rich_visibility, &error)) << error;
  EXPECT_TRUE(error.empty());
  ASSERT_EQ(rich_visibility.size(), first_rich_visibility.size());
  for (size_t i = 0; i < rich_visibility.size(); ++i) {
    EXPECT_EQ(rich_visibility[i].position, first_rich_visibility[i].position);
    EXPECT_TRUE(rich_visibility[i].support == first_rich_visibility[i].support);
    EXPECT_EQ(exactPoint(rich_visibility[i]), exactPoint(first_rich_visibility[i]));
  }

  std::vector<std::pair<double, double>> legacy_visibility;
  std::vector<std::pair<int, int>> legacy_topology;
  ASSERT_TRUE(poly.getVisibilityRegion(5, 5, legacy_visibility, legacy_topology, &error)) << error;
  ASSERT_EQ(legacy_visibility.size(), rich_visibility.size());
  ASSERT_EQ(legacy_topology.size(), rich_visibility.size());
  for (size_t i = 0; i < rich_visibility.size(); ++i) {
    EXPECT_EQ(legacy_visibility[i], rich_visibility[i].position);
    if (const auto* vertex = std::get_if<ObstacleVertexId>(&rich_visibility[i].support)) {
      EXPECT_EQ(legacy_topology[i], std::make_pair(vertex->obstacle, vertex->vertex));
    } else {
      EXPECT_EQ(legacy_topology[i], std::make_pair(-1, -1));
    }
  }
}

TEST(Polymap, ContinuousRootVisibilityIsSourceSpecificAndUncached) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);
  ASSERT_TRUE(poly.hasSolution()) << PolymapTestPeer::constructionError(poly);
  const size_t cache_before = PolymapTestPeer::visibilityCacheSize(poly);

  VisibilityRegion first;
  VisibilityRegion second;
  std::string error;
  ASSERT_TRUE(poly.getRootVisibilityRegion(Point2d{5.25, 5.25}, first, &error)) << error;
  ASSERT_TRUE(poly.getRootVisibilityRegion(Point2d{5.75, 5.75}, second, &error)) << error;
  ASSERT_GE(first.size(), 2u);
  ASSERT_GE(second.size(), 2u);
  EXPECT_TRUE(poly.validateVisibilityRegion(first, &error)) << error;
  EXPECT_TRUE(poly.validateVisibilityRegion(second, &error)) << error;
  EXPECT_EQ(PolymapTestPeer::visibilityCacheSize(poly), cache_before);

  bool differs = first.size() != second.size();
  if (!differs) {
    for (size_t i = 0; i < first.size(); ++i) {
      if (first[i].position != second[i].position || !(first[i].support == second[i].support)) {
        differs = true;
        break;
      }
    }
  }
  EXPECT_TRUE(differs);
}

TEST(Polymap, CooperativeVisibilityStopDoesNotCommitOutputOrCache) {
  const auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);
  ASSERT_TRUE(poly.isCDTReady()) << PolymapTestPeer::constructionError(poly);
  ASSERT_EQ(PolymapTestPeer::visibilityCacheSize(poly), 0u);

  VisibilityRegion baseline;
  std::string error;
  size_t full_poll_count = 0;
  const StopToken count_only([&full_poll_count]() {
    ++full_poll_count;
    return false;
  });
  ASSERT_EQ(poly.getVisibilityRegion(5, 5, baseline, count_only, &error), OperationStatus::success)
    << error;
  ASSERT_FALSE(baseline.empty());
  ASSERT_GT(full_poll_count, 0u);
  ASSERT_EQ(PolymapTestPeer::visibilityCacheSize(poly), 1u);
  PolymapTestPeer::clearVisibilityCache(poly);

  VisibilityRegion visibility = {BoundaryEndpoint{{-100.0, -100.0}, std::monostate{}}};
  error = "stale error";
  size_t stopped_poll_count = 0;
  const StopToken stop_at_last_poll(
    [&stopped_poll_count, full_poll_count]() { return ++stopped_poll_count >= full_poll_count; });

  EXPECT_EQ(poly.getVisibilityRegion(5, 5, visibility, stop_at_last_poll, &error),
            OperationStatus::stopped);
  EXPECT_EQ(stopped_poll_count, full_poll_count);
  EXPECT_TRUE(visibility.empty());
  EXPECT_TRUE(error.empty());
  EXPECT_EQ(PolymapTestPeer::visibilityCacheSize(poly), 0u);

  EXPECT_TRUE(poly.getVisibilityRegion(5, 5, visibility, &error)) << error;
  ASSERT_EQ(visibility.size(), baseline.size());
  for (size_t index = 0; index < baseline.size(); ++index) {
    EXPECT_EQ(visibility[index].position, baseline[index].position);
    EXPECT_EQ(exactPoint(visibility[index]), exactPoint(baseline[index]));
    EXPECT_EQ(visibility[index].support, baseline[index].support);
  }
  EXPECT_EQ(PolymapTestPeer::visibilityCacheSize(poly), 1u);
}

TEST(Polymap, StoppedCachedVisibilityValidationKeepsValidCache) {
  const auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);
  ASSERT_TRUE(poly.isCDTReady()) << PolymapTestPeer::constructionError(poly);

  VisibilityRegion expected;
  std::string error;
  ASSERT_TRUE(poly.getVisibilityRegion(5, 5, expected, &error)) << error;
  ASSERT_FALSE(expected.empty());
  ASSERT_EQ(PolymapTestPeer::visibilityCacheSize(poly), 1u);

  VisibilityRegion cached_baseline;
  size_t full_poll_count = 0;
  const StopToken count_only([&full_poll_count]() {
    ++full_poll_count;
    return false;
  });
  ASSERT_EQ(poly.getVisibilityRegion(5, 5, cached_baseline, count_only, &error),
            OperationStatus::success)
    << error;
  ASSERT_GT(full_poll_count, 0u);

  VisibilityRegion stopped_output = {BoundaryEndpoint{{-100.0, -100.0}, std::monostate{}}};
  error = "stale error";
  size_t stopped_poll_count = 0;
  const StopToken stop_at_last_poll(
    [&stopped_poll_count, full_poll_count]() { return ++stopped_poll_count >= full_poll_count; });

  EXPECT_EQ(poly.getVisibilityRegion(5, 5, stopped_output, stop_at_last_poll, &error),
            OperationStatus::stopped);
  EXPECT_EQ(stopped_poll_count, full_poll_count);
  EXPECT_TRUE(stopped_output.empty());
  EXPECT_TRUE(error.empty());
  EXPECT_EQ(PolymapTestPeer::visibilityCacheSize(poly), 1u);

  VisibilityRegion recovered;
  EXPECT_TRUE(poly.getVisibilityRegion(5, 5, recovered, &error)) << error;
  ASSERT_EQ(recovered.size(), cached_baseline.size());
  for (size_t index = 0; index < cached_baseline.size(); ++index) {
    EXPECT_EQ(recovered[index].position, cached_baseline[index].position);
    EXPECT_EQ(exactPoint(recovered[index]), exactPoint(cached_baseline[index]));
    EXPECT_EQ(recovered[index].support, cached_baseline[index].support);
  }
  EXPECT_EQ(PolymapTestPeer::visibilityCacheSize(poly), 1u);
}

TEST(Polymap, ObstacleVertexVisibilityUsesOpenSectorAnchors) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);
  ASSERT_TRUE(poly.hasSolution());

  constexpr int source_x = 10;
  constexpr int source_y = 10;
  const auto source_topology = poly.locateVertex(source_x, source_y);
  ASSERT_TRUE(poly.isValidTopology(source_topology));
  const auto previous = poly.getPrevObs(source_topology);
  const auto next = poly.getNextObs(source_topology);
  ASSERT_TRUE(previous.has_value());
  ASSERT_TRUE(next.has_value());

  VisibilityRegion visibility;
  std::string error;
  ASSERT_TRUE(poly.getVisibilityRegion(source_x, source_y, visibility, &error)) << error;
  ASSERT_GE(visibility.size(), 2u);
  EXPECT_EQ(exactPoint(visibility.front()),
            exact_geometry::Point(previous->first, previous->second));
  EXPECT_EQ(exactPoint(visibility.back()), exact_geometry::Point(next->first, next->second));

  const exact_geometry::Point source(source_x, source_y);
  for (const auto& endpoint : visibility) EXPECT_NE(exactPoint(endpoint), source);
}

TEST(Polymap, CachedVisibilityIsRevalidatedWithSourceGeometry) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);
  ASSERT_TRUE(poly.hasSolution());

  VisibilityRegion visibility;
  std::string error;
  ASSERT_TRUE(poly.getVisibilityRegion(5, 5, visibility, &error)) << error;
  ASSERT_GE(visibility.size(), 3u);
  ASSERT_TRUE(PolymapTestPeer::reverseCachedVisibility(poly, 5, 5));

  visibility = {BoundaryEndpoint{{-100.0, -100.0}, std::monostate{}}};
  error = "stale error";
  EXPECT_FALSE(poly.getVisibilityRegion(5, 5, visibility, &error));
  EXPECT_TRUE(visibility.empty());
  EXPECT_NE(error.find("Cached"), std::string::npos) << error;
  EXPECT_NE(error.find("ray order"), std::string::npos) << error;
  EXPECT_EQ(PolymapTestPeer::visibilityCacheSize(poly), 0u);

  error.clear();
  EXPECT_TRUE(poly.getVisibilityRegion(5, 5, visibility, &error)) << error;
  EXPECT_FALSE(visibility.empty());
}

TEST(Polymap, VisibilityFailureClearsOutputsOnEveryCall) {
  GridMap map;
  map.width = 20;
  map.height = 20;
  map.resolution = 1.0f;
  map.data.resize(400, 0);

  Polymap poly(map, 2, 2, 17, 17);
  ASSERT_TRUE(poly.hasSolution());
  ASSERT_TRUE(poly.obstacles().empty());

  std::vector<std::pair<double, double>> visibility = {{1.0, 2.0}};
  std::vector<std::pair<int, int>> topology = {{3, 4}};
  std::string error;
  EXPECT_FALSE(poly.getVisibilityRegion(2, 2, visibility, topology, &error));
  EXPECT_TRUE(visibility.empty());
  EXPECT_TRUE(topology.empty());
  EXPECT_FALSE(error.empty());

  visibility = {{5.0, 6.0}};
  topology = {{7, 8}, {9, 10}};
  error.clear();
  EXPECT_FALSE(poly.getVisibilityRegion(2, 2, visibility, topology, &error));
  EXPECT_TRUE(visibility.empty());
  EXPECT_TRUE(topology.empty());
  EXPECT_FALSE(error.empty());
}

TEST(Polymap, NonFiniteRichEndpointIsRejectedWithoutExactConstruction) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);
  ASSERT_TRUE(poly.hasSolution());

  const BoundaryEndpoint endpoint{{std::numeric_limits<double>::quiet_NaN(), 5.0},
                                  std::monostate{}};
  EXPECT_FALSE(endpoint.exact_position.has_value());

  std::string error;
  EXPECT_FALSE(poly.validateBoundaryEndpoint(endpoint, &error));
  EXPECT_NE(error.find("not finite"), std::string::npos);
}

TEST(Polymap, TopologyAccessorsRejectInvalidIndices) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);

  ASSERT_TRUE(poly.hasSolution());
  ASSERT_FALSE(poly.obstacles().empty());
  ASSERT_GE(poly.obstacles()[0].ordered_vertices_.size(), 2u);

  const int vertex_count = static_cast<int>(poly.obstacles()[0].ordered_vertices_.size());
  EXPECT_TRUE(poly.isValidTopology({0, 0}));
  EXPECT_FALSE(poly.isValidTopology({-1, -1}));
  EXPECT_FALSE(poly.isValidTopology({-1, 0}));
  EXPECT_FALSE(poly.isValidTopology({0, -1}));
  EXPECT_FALSE(poly.isValidTopology({static_cast<int>(poly.obstacles().size()), 0}));
  EXPECT_FALSE(poly.isValidTopology({0, vertex_count}));

  EXPECT_FALSE(poly.areConsecutive({-1, -1}, {-1, -1}));
  EXPECT_FALSE(poly.areConsecutive({0, -1}, {0, 0}));
  EXPECT_FALSE(poly.areConsecutive({0, 0}, {0, vertex_count}));
  EXPECT_FALSE(poly.getPrevObs({0, -1}).has_value());
  EXPECT_FALSE(poly.getNextObs({0, vertex_count}).has_value());

  auto prev = poly.getPrevObs({0, 0});
  auto next = poly.getNextObs({0, 0});
  ASSERT_TRUE(prev.has_value());
  ASSERT_TRUE(next.has_value());
  EXPECT_EQ(*prev, poly.obstacles()[0].ordered_vertices_.back());
  EXPECT_EQ(*next, poly.obstacles()[0].ordered_vertices_[1]);
  EXPECT_TRUE(poly.areConsecutive({0, 0}, {0, 1}));
  EXPECT_TRUE(poly.areConsecutive({0, vertex_count - 1}, {0, 0}));

  poly.obstacles().emplace_back(Obs());
  const int empty_obstacle = static_cast<int>(poly.obstacles().size()) - 1;
  EXPECT_FALSE(poly.isValidTopology({empty_obstacle, 0}));
  EXPECT_FALSE(poly.areConsecutive({empty_obstacle, 0}, {empty_obstacle, 0}));
  EXPECT_FALSE(poly.getPrevObs({empty_obstacle, 0}).has_value());
  EXPECT_FALSE(poly.getNextObs({empty_obstacle, 0}).has_value());
}

TEST(Polymap, RichBoundarySupportIsDirectionalAndCannotBeDropped) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);
  ASSERT_TRUE(poly.hasSolution());

  DirectedObstacleEdge edge;
  std::pair<int, int> from;
  std::pair<int, int> to;
  bool found_long_edge = false;
  for (size_t obstacle_index = 0; obstacle_index < poly.obstacles().size() && !found_long_edge;
       ++obstacle_index) {
    const auto& vertices = poly.obstacles()[obstacle_index].ordered_vertices_;
    for (size_t vertex_index = 0; vertex_index < vertices.size(); ++vertex_index) {
      const size_t next_index = (vertex_index + 1) % vertices.size();
      const double length = std::hypot(
        static_cast<double>(vertices[next_index].first - vertices[vertex_index].first),
        static_cast<double>(vertices[next_index].second - vertices[vertex_index].second));
      if (length <= 2.0)
        continue;
      edge = DirectedObstacleEdge{
        ObstacleVertexId{static_cast<int>(obstacle_index), static_cast<int>(vertex_index)},
        ObstacleVertexId{static_cast<int>(obstacle_index), static_cast<int>(next_index)}};
      from = vertices[vertex_index];
      to = vertices[next_index];
      found_long_edge = true;
      break;
    }
  }
  ASSERT_TRUE(found_long_edge);

  const auto interpolate = [&](double parameter) {
    return Point2d{from.first + parameter * (to.first - from.first),
                   from.second + parameter * (to.second - from.second)};
  };
  const BoundaryEndpoint exact_from{
    {static_cast<double>(from.first), static_cast<double>(from.second)}, edge.from};
  const BoundaryEndpoint exact_to{{static_cast<double>(to.first), static_cast<double>(to.second)},
                                  edge.to};
  const BoundaryEndpoint first_interior{interpolate(0.25), edge};
  const BoundaryEndpoint second_interior{interpolate(0.75), edge};
  const BoundaryEndpoint missing_support{interpolate(0.5), std::monostate{}};

  std::string error;
  EXPECT_FALSE(poly.validateBoundaryEndpoint(missing_support, &error));
  EXPECT_NE(error.find("supporting-edge"), std::string::npos);
  EXPECT_TRUE(poly.boundarySupportsConsecutive(exact_from, exact_to));
  EXPECT_FALSE(poly.boundarySupportsConsecutive(exact_to, exact_from));
  EXPECT_TRUE(poly.boundarySupportsConsecutive(exact_from, first_interior));
  EXPECT_TRUE(poly.boundarySupportsConsecutive(first_interior, second_interior));
  EXPECT_FALSE(poly.boundarySupportsConsecutive(second_interior, first_interior));
  EXPECT_TRUE(poly.boundarySupportsConsecutive(second_interior, exact_to));
  EXPECT_FALSE(
    poly.boundarySupportsConsecutive(BoundaryEndpoint{{5.5, 5.5}, std::monostate{}}, exact_from));
}

TEST(Polymap, LocateVertexRejectsInvalidCoordinates) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);

  EXPECT_EQ(poly.locateVertex(-1, 0), std::make_pair(-1, -1));
  EXPECT_EQ(poly.locateVertex(0, -1), std::make_pair(-1, -1));
  EXPECT_EQ(poly.locateVertex(static_cast<int>(map.width), 0), std::make_pair(-1, -1));
  EXPECT_EQ(poly.locateVertex(0, static_cast<int>(map.height)), std::make_pair(-1, -1));
  EXPECT_EQ(poly.locateVertex(1.25, 1.0), std::make_pair(-1, -1));
  EXPECT_EQ(poly.locateVertex(std::numeric_limits<double>::quiet_NaN(), 1.0),
            std::make_pair(-1, -1));
  EXPECT_EQ(poly.locateVertex(std::numeric_limits<double>::infinity(), 1.0),
            std::make_pair(-1, -1));
}

TEST(Polymap, BoundedCDTEdgeExportDoesNotCopyTheWholeTriangulation) {
  auto map = makeMapWithObstacle();
  Polymap poly(map, 5, 5, 25, 25);
  ASSERT_TRUE(poly.hasSolution());
  ASSERT_TRUE(poly.isCDTReady());

  const auto all_edges = poly.getCDTEdges();
  const auto limited_edges = poly.getCDTEdges(3);
  EXPECT_LE(limited_edges.size(), 3u);
  EXPECT_LE(limited_edges.size(), all_edges.size());
  EXPECT_TRUE(poly.getCDTEdges(0).empty());
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
