#pragma once

#include <cstddef>
#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <raystar/coordinate_utils.h>
#include <raystar/planning_limits.h>
#include <raystar/polymap.h>

namespace raystar
{

class Candidate
{
public:
  Candidate(int node_index, int child_index, double cost)
    : Nindex_(node_index), Cindex_(child_index), Fcost_(cost) {}

  Candidate(const Candidate&) = default;

  [[nodiscard]] int nodeIndex() const noexcept { return Nindex_; }
  [[nodiscard]] int childIndex() const noexcept { return Cindex_; }
  [[nodiscard]] double fCost() const noexcept { return Fcost_; }

#ifdef RAYSTAR_TESTING
public:  // Legacy white-box tests; production builds keep these fields private.
#else
private:
#endif
  friend class RaystarCore;

  int Nindex_;
  int Cindex_;
  double Fcost_;
};

class Child
{
public:
  Child(int nindex, int cindex, int cx, int cy, bool is_left_gap)
    : Nindex_(nindex), Cindex_(cindex), start_angle_(0), end_angle_(0),
      c_({cx, cy}), o_({0.0, 0.0}),
      c_endpoint_({{static_cast<double>(cx), static_cast<double>(cy)}, std::monostate{}}),
      o_endpoint_(),
      c_obs_index_(-1), c_ver_index_(-1), o_obs_index_(-1), o_ver_index_(-1),
      is_a_left_gap_(is_left_gap), c_gcost_(0), c_hcost_(0) {}

  Child(int nindex, int cindex, const BoundaryEndpoint& c_endpoint,
    const BoundaryEndpoint& o_endpoint, bool is_left_gap);

  // Construct a fully initialized visibility-gap child without exposing
  // mutable setters for its search metadata.
  Child(int nindex, int cindex, const BoundaryEndpoint& c_endpoint,
    const BoundaryEndpoint& o_endpoint, bool is_left_gap,
    double start_angle, double end_angle, double g_cost);

  [[nodiscard]] int nodeIndex() const noexcept { return Nindex_; }
  [[nodiscard]] int childIndex() const noexcept { return Cindex_; }
  [[nodiscard]] double startAngle() const noexcept { return start_angle_; }
  [[nodiscard]] double endAngle() const noexcept { return end_angle_; }
  [[nodiscard]] const std::pair<int, int>& coordinate() const noexcept { return c_; }
  [[nodiscard]] const Point2d& oppositeCoordinate() const noexcept { return o_; }
  [[nodiscard]] const BoundaryEndpoint& endpoint() const noexcept
  {
    return c_endpoint_;
  }
  [[nodiscard]] const BoundaryEndpoint& oppositeEndpoint() const noexcept
  {
    return o_endpoint_;
  }
  [[nodiscard]] int obstacleIndex() const noexcept { return c_obs_index_; }
  [[nodiscard]] int vertexIndex() const noexcept { return c_ver_index_; }
  [[nodiscard]] int oppositeObstacleIndex() const noexcept { return o_obs_index_; }
  [[nodiscard]] int oppositeVertexIndex() const noexcept { return o_ver_index_; }
  [[nodiscard]] bool isLeftGap() const noexcept { return is_a_left_gap_; }
  [[nodiscard]] double gCost() const noexcept { return c_gcost_; }
  [[nodiscard]] double hCost() const noexcept { return c_hcost_; }

#ifdef RAYSTAR_TESTING
public:  // Legacy white-box tests; production builds keep these fields private.
#else
private:
#endif
  friend class Node;
  friend class RaystarCore;
  friend class NodeTestPeer;

  int Nindex_;
  int Cindex_;
  double start_angle_;
  double end_angle_;

  std::pair<int, int> c_;
  Point2d o_;
  BoundaryEndpoint c_endpoint_;
  BoundaryEndpoint o_endpoint_;

  int c_obs_index_;
  int c_ver_index_;
  int o_obs_index_;
  int o_ver_index_;

  bool is_a_left_gap_;

  double c_gcost_;
  double c_hcost_;
};

class Node
{
public:
  Node(const Polymap* pMap, int Nindex, double start_x, double start_y,
    double Gcost, double Hcost,
    const VisibilityRegion& visibility_region);

  Node(const Polymap* pMap, int Nindex, double start_x, double start_y,
    double Gcost, double Hcost,
    const VisibilityRegion& visibility_region,
    const StopToken& stop_token);

  Node(const Polymap* pMap, int Nindex, double seed_x, double seed_y,
    double Gcost, double Hcost, int parent_index,
    const VisibilityRegion& visibility_region);

  Node(const Polymap* pMap, int Nindex, double seed_x, double seed_y,
    double Gcost, double Hcost, int parent_index,
    const VisibilityRegion& visibility_region,
    const StopToken& stop_token);

  Node(const Polymap* pMap, int Nindex, double start_x, double start_y,
    double Gcost, double Hcost,
    const std::vector<std::pair<double, double>>& visibility_region,
    const std::vector<std::pair<int, int>>& topo_V);

  Node(const Polymap* pMap, int Nindex, double seed_x, double seed_y,
    double Gcost, double Hcost, int parent_index,
    const std::vector<std::pair<double, double>>& visibility_region,
    const std::vector<std::pair<int, int>>& topo_V);

  void setFullVisibilityRegion(const VisibilityRegion& visibility_region);
  [[nodiscard]] OperationStatus setFullVisibilityRegion(
    const VisibilityRegion& visibility_region,
    const StopToken& stop_token);

  [[nodiscard]] bool generateChild(const Polymap* pMap);
  [[nodiscard]] bool generateChild(const Polymap* pMap, std::string* error);
  [[nodiscard]] OperationStatus generateChild(
    const Polymap* pMap, const StopToken& stop_token,
    std::string* error = nullptr);

  [[nodiscard]] int index() const noexcept { return Nindex_; }
  [[nodiscard]] const std::pair<int, int>& seed() const noexcept { return seed_; }
  [[nodiscard]] bool hasValidSeed() const noexcept { return seed_is_valid_; }
  [[nodiscard]] bool isContinuousRoot() const noexcept { return is_continuous_root_; }
  [[nodiscard]] double startAngle() const noexcept { return start_angle_; }
  [[nodiscard]] double endAngle() const noexcept { return end_angle_; }
  [[nodiscard]] int parentIndex() const noexcept { return parent_index_; }
  [[nodiscard]] const Point2d& startBoundary() const noexcept { return start_o_; }
  [[nodiscard]] const Point2d& endBoundary() const noexcept { return end_o_; }
  [[nodiscard]] bool isChildLeftGap() const noexcept { return as_a_child_left_gap_; }
  [[nodiscard]] double gCost() const noexcept { return Gcost_; }
  [[nodiscard]] double hCost() const noexcept { return Hcost_; }
  [[nodiscard]] const std::vector<Child>& children() const noexcept { return C_; }
  [[nodiscard]] const std::vector<Point2d>& visibility() const noexcept { return V_; }
  [[nodiscard]] const std::vector<std::pair<int, int>>& visibilityTopology() const noexcept
  {
    return topo_V_;
  }
  [[nodiscard]] const std::vector<Point2d>& fullVisibility() const noexcept
  {
    return full_V_;
  }
  [[nodiscard]] const VisibilityRegion& visibilityRegion() const noexcept
  {
    return visibility_region_;
  }
  [[nodiscard]] const VisibilityRegion& fullVisibilityRegion() const noexcept
  {
    return full_visibility_region_;
  }
  [[nodiscard]] bool hasValidVisibilityRegion() const noexcept
  {
    return visibility_region_valid_;
  }
  [[nodiscard]] const std::string& visibilityRegionError() const noexcept
  {
    return visibility_region_error_;
  }
  [[nodiscard]] const std::vector<std::pair<int, int>>& localShortestPath() const noexcept
  {
    return local_shortest_path_;
  }
  [[nodiscard]] const std::vector<int>& pathNodeIndices() const noexcept
  {
    return path_node_index_;
  }

#ifdef RAYSTAR_TESTING
public:  // Legacy white-box tests; production builds keep these fields private.
#else
private:
#endif
  friend class RaystarCore;
  friend class NodeTestPeer;

  int Nindex_;
  std::pair<int, int> seed_;
  bool seed_is_valid_;
  // True only for the compatibility/debug root holder created by the
  // continuous-endpoint planner. Its integer seed is not root geometry.
  bool is_continuous_root_ = false;
  double start_angle_;
  double end_angle_;
  int parent_index_;

  Point2d start_o_;
  Point2d end_o_;
  bool as_a_child_left_gap_;

  double Gcost_;
  double Hcost_;

  std::vector<Child> C_;

  std::vector<Point2d> V_;
  std::vector<std::pair<int, int>> topo_V_;
  std::vector<Point2d> full_V_;
  VisibilityRegion visibility_region_;
  VisibilityRegion full_visibility_region_;
  bool visibility_region_valid_;
  std::string visibility_region_error_;

  std::vector<std::pair<int, int>> local_shortest_path_;
  std::vector<int> path_node_index_;
};

/// A planning endpoint has two intentionally different representations.
/// `cell_` is used only by the occupancy grid/flood fill; `position_` is the
/// continuous grid coordinate used by visibility, distance and path geometry.
struct PlanEndpoint
{
  std::pair<int, int> cell_ = {0, 0};
  Point2d position_ = {0.0, 0.0};

  PlanEndpoint() = default;
  PlanEndpoint(std::pair<int, int> cell, Point2d position)
    : cell_(cell), position_(position) {}
  PlanEndpoint(int cell_x, int cell_y, double x, double y)
    : cell_({cell_x, cell_y}), position_({x, y}) {}
};

struct PathSolution
{
  // Endpoints remain continuous grid coordinates.  Only turning_points_ are
  // stored as integer obstacle vertices; projectedPath() is a temporary view
  // for ROS/debug consumers.
  Point2d start_;
  std::vector<std::pair<int, int>> turning_points_;
  Point2d goal_;
  // Source-compatible view retained for callers of the original API.  It is
  // an intentionally lossy cell projection: continuous endpoints are mapped
  // with floor(), while turning points are copied exactly.  New code should
  // use start_/turning_points_/goal_ or projectedPath() instead.
  std::vector<std::pair<int, int>> path_;
  double path_cost_ = 0.0;
  std::vector<int> path_node_index_;

  PathSolution() = default;
  // Compatibility constructor used by the original integer-only API.
  PathSolution(const std::vector<std::pair<int, int>>& path,
    double path_cost, const std::vector<int>& path_node_index)
    : path_(path), path_cost_(path_cost)
  {
    if (!path.empty())
    {
      start_ = {static_cast<double>(path.front().first),
        static_cast<double>(path.front().second)};
      goal_ = start_;
      if (path.size() > 1)
      {
        goal_ = {static_cast<double>(path.back().first),
          static_cast<double>(path.back().second)};
      }
      if (path.size() > 2)
      {
        turning_points_.assign(path.begin() + 1, path.end() - 1);
      }
    }
    path_node_index_.assign(path_node_index.begin(), path_node_index.end());
  }

  PathSolution(const Point2d& start,
    const std::vector<std::pair<int, int>>& turning_points,
    const Point2d& goal, double path_cost,
    const std::vector<int>& path_node_index)
    : start_(start), turning_points_(turning_points), goal_(goal),
      path_cost_(path_cost)
  {
    const auto legacy_coordinate = [](double coordinate) {
      if (!std::isfinite(coordinate) ||
          coordinate < static_cast<double>(std::numeric_limits<int>::lowest()) ||
          coordinate > static_cast<double>(std::numeric_limits<int>::max()))
      {
        return 0;
      }
      return static_cast<int>(std::floor(coordinate));
    };
    path_.reserve(turning_points_.size() + 2);
    path_.emplace_back(
      legacy_coordinate(start_.first), legacy_coordinate(start_.second));
    path_.insert(path_.end(), turning_points_.begin(), turning_points_.end());
    path_.emplace_back(
      legacy_coordinate(goal_.first), legacy_coordinate(goal_.second));
    path_node_index_.assign(path_node_index.begin(), path_node_index.end());
  }

  [[nodiscard]] std::vector<Point2d> projectedPath() const
  {
    std::vector<Point2d> result;
    result.reserve(turning_points_.size() + 2);
    result.emplace_back(start_);
    for (const auto& point : turning_points_)
    {
      result.emplace_back(
        static_cast<double>(point.first), static_cast<double>(point.second));
    }
    result.emplace_back(goal_);
    return result;
  }
};

enum class PlanningLimitReached
{
  none,
  max_nodes,
  timeout,
  cancelled,
  max_path_points
};

// Describes how Core reached its terminal result without requiring ROS or
// direct callers to parse the diagnostic message. `complete` means either K
// paths were found or the search frontier was exhausted with at least one
// path; `no_path` is the corresponding exhausted zero-path result.
enum class PlanningOutcome
{
  invalid_request,
  complete,
  no_path,
  limit_reached,
  failed
};

struct PlanResult
{
  bool success = false;
  std::string message;
  std::vector<PathSolution> path_solutions;
  PlanningOutcome outcome = PlanningOutcome::failed;
  PlanningLimitReached limit_reached = PlanningLimitReached::none;
  size_t expanded_nodes = 0;
  double map_time_ms = 0.0;
  double plan_time_ms = 0.0;
  std::shared_ptr<const Polymap> polymap;  // owns the Polymap; safe to use after plan() returns
};

class RaystarCore
{
public:
  RaystarCore() = default;

  PlanResult plan(const GridMap& grid_map,
    int start_x, int start_y, int goal_x, int goal_y,
    int K, bool allow_self_crossing,
    const PlanningLimits& limits = PlanningLimits{});

  PlanResult plan(const GridMap& grid_map,
    const Point2d& start, const Point2d& goal,
    int K, bool allow_self_crossing,
    const PlanningLimits& limits = PlanningLimits{});

  PlanResult plan(const GridMap& grid_map,
    const PlanEndpoint& start, const PlanEndpoint& goal,
    int K, bool allow_self_crossing,
    const PlanningLimits& limits = PlanningLimits{});

  const std::vector<Node>& getNodes() const { return N_; }

  // Release search-tree objects (including their vector capacities) between
  // requests.  The ROS node uses this before admitting a new map so a large
  // previous request cannot remain resident while the next request is being
  // copied and converted.
  void resetSearchState() noexcept
  {
    std::vector<Node>().swap(N_);
  }

private:
  void outlineMap(std::vector<uint8_t>& costarr, int nx, int ny);
  [[nodiscard]] OperationStatus outlineMap(
    std::vector<uint8_t>& costarr, int nx, int ny,
    const StopToken& stop_token);
  [[nodiscard]] OperationStatus getScopedVisibilityRegion(Polymap& theMap,
    Candidate& the_child,
    VisibilityRegion& visibility_region,
    VisibilityRegion& full_visreg,
    const StopToken& stop_token, std::string& error);

  std::vector<Node> N_;
};

}  // namespace raystar
