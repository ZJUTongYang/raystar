#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <optional>
#include <limits>
#include <string>
#include <utility>
#include <variant>
#include <raystar/coordinate_utils.h>
#include <raystar/exact_geometry.h>
#include <raystar/planning_limits.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <raystar/cooperative_stop.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

namespace raystar
{

using Point2d = std::pair<double, double>;

struct ObstacleVertexId
{
  int obstacle = -1;
  int vertex = -1;

  bool operator==(const ObstacleVertexId& other) const
  {
    return obstacle == other.obstacle && vertex == other.vertex;
  }
};

struct DirectedObstacleEdge
{
  ObstacleVertexId from;
  ObstacleVertexId to;

  bool operator==(const DirectedObstacleEdge& other) const
  {
    return from == other.from && to == other.to;
  }
};

using BoundarySupport = std::variant<
  std::monostate,
  ObstacleVertexId,
  DirectedObstacleEdge>;

struct BoundaryEndpoint
{
  Point2d position = {0.0, 0.0};
  BoundarySupport support;
  // Internal geometry keeps the EPECK construction. `position` remains the
  // public/ROS double projection; legacy callers may omit exact_position, in
  // which case the supplied doubles are treated as exact binary rationals.
  std::optional<exact_geometry::Point> exact_position;

  BoundaryEndpoint() = default;

  BoundaryEndpoint(Point2d projected_position, BoundarySupport boundary_support)
    : position(projected_position), support(std::move(boundary_support))
  {
    // Do not pass NaN/Inf into CGAL. Validation will report the endpoint as a
    // normal input error before exactPoint() is used by planning code.
    if (std::isfinite(projected_position.first) &&
        std::isfinite(projected_position.second))
    {
      exact_position.emplace(
        projected_position.first, projected_position.second);
    }
  }
};

inline exact_geometry::Point exactPoint(const BoundaryEndpoint& endpoint)
{
  if (endpoint.exact_position)
    return *endpoint.exact_position;
  return exact_geometry::Point(endpoint.position.first, endpoint.position.second);
}

using VisibilityRegion = std::vector<BoundaryEndpoint>;

namespace constrained_delaunay_triangulation
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, CGAL::No_constraint_intersection_tag> CDT;
  typedef CDT::Point Point;

  struct BungiuEdge
  {
    BoundaryEndpoint prev;
    BoundaryEndpoint next;
    double limit_prev = 0.0;
    double limit_next = 0.0;
    BoundaryEndpoint limit_prev_anchor;
    BoundaryEndpoint limit_next_anchor;
    bool is_bd = false;
    std::optional<DirectedObstacleEdge> supporting_segment;
  };
}

class Obs
{
public:
  std::vector<std::pair<int, int>> ordered_vertices_;
};

class PolymapTestPeer;
class Polymap;

enum class PolymapCreateStatus
{
  ready,
  no_path,
  stopped,
  failure
};

struct PolymapCreateResult;

class Polymap
{
public:
  Polymap(const Polymap&) = delete;
  Polymap& operator=(const Polymap&) = delete;
  Polymap(Polymap&&) = default;
  Polymap& operator=(Polymap&&) = default;
  ~Polymap() = default;

  // Controlled construction never exposes a partially initialized Polymap.
  // `ready` owns the value; all other statuses carry only a diagnostic.
  [[nodiscard]] static PolymapCreateResult create(
    const GridMap& grid_map, int start_x, int start_y,
    int goal_x, int goal_y,
    const PlanningLimits& limits = PlanningLimits{});
  [[nodiscard]] static PolymapCreateResult create(
    const GridMap& grid_map, int start_x, int start_y,
    int goal_x, int goal_y, const StopToken& stop_token,
    const PlanningLimits& limits = PlanningLimits{});
  [[nodiscard]] static PolymapCreateResult create(
    const GridMap& grid_map, int start_x, int start_y,
    int goal_x, int goal_y, const Point2d& start_position,
    const Point2d& goal_position,
    const PlanningLimits& limits = PlanningLimits{});
  [[nodiscard]] static PolymapCreateResult create(
    const GridMap& grid_map, int start_x, int start_y,
    int goal_x, int goal_y, const Point2d& start_position,
    const Point2d& goal_position, const StopToken& stop_token,
    const PlanningLimits& limits = PlanningLimits{});

  [[nodiscard]] bool isCDTReady() const noexcept { return cdt_ready_; }
  [[nodiscard]] bool hasSolution() const noexcept { return solution_exist_; }
  [[nodiscard]] int width() const noexcept { return xsize_; }
  [[nodiscard]] int height() const noexcept { return ysize_; }
  [[nodiscard]] const std::vector<uint8_t>& occupancyData() const noexcept
  {
    return data_;
  }
  [[nodiscard]] const std::vector<Obs>& obstacles() const noexcept
  {
    return obs_;
  }
#ifdef RAYSTAR_TESTING
  // Legacy white-box diagnostics are deliberately unavailable in production;
  // Polymap::create() returns status/error without exposing a partially built
  // object's construction state as public API.
  [[nodiscard]] bool constructionStopped() const noexcept
  {
    return construction_stopped_;
  }
  [[nodiscard]] const std::string& constructionError() const noexcept
  {
    return construction_error_;
  }
  [[nodiscard]] std::vector<Obs>& obstacles() noexcept
  {
    return obs_;
  }
#endif

  [[nodiscard]] bool getVisibilityRegion(int start_x, int start_y,
    VisibilityRegion& visibility_region);

  [[nodiscard]] bool getVisibilityRegion(int start_x, int start_y,
    VisibilityRegion& visibility_region,
    std::string* error);

  [[nodiscard]] OperationStatus getVisibilityRegion(
    int start_x, int start_y, VisibilityRegion& visibility_region,
    const StopToken& stop_token, std::string* error = nullptr);

  // Continuous free-space sources are used only for the root expansion.  This
  // entry point deliberately bypasses the integer obstacle-vertex cache: two
  // positions in the same grid cell generally have different visibility.
  [[nodiscard]] bool getRootVisibilityRegion(
    const Point2d& source, VisibilityRegion& visibility_region,
    std::string* error = nullptr);
  [[nodiscard]] OperationStatus getRootVisibilityRegion(
    const Point2d& source, VisibilityRegion& visibility_region,
    const StopToken& stop_token, std::string* error = nullptr);

  // Exact, epsilon-free validation against the final simplified free-space
  // boundary.  A valid point must be strictly inside the reachable outer
  // contour and strictly outside every obstacle hole.
  [[nodiscard]] bool validateFreeSpaceInterior(
    const Point2d& point, std::string* error = nullptr) const;
  [[nodiscard]] OperationStatus validateFreeSpaceInterior(
    const Point2d& point, const StopToken& stop_token,
    std::string* error = nullptr) const;

  // Legacy, intentionally lossy compatibility API. Only exact obstacle
  // vertices receive a topology index; edge-interior intersections are
  // reported as (-1, -1). Planning code must use the rich overloads above.
  [[nodiscard]] bool getVisibilityRegion(int start_x, int start_y,
    std::vector<std::pair<double, double>>& visibility_region,
    std::vector<std::pair<int, int>>& topo_V);

  [[nodiscard]] bool getVisibilityRegion(int start_x, int start_y,
    std::vector<std::pair<double, double>>& visibility_region,
    std::vector<std::pair<int, int>>& topo_V,
    std::string* error);

  std::pair<int, int> locateVertex(int x, int y) const;
  std::pair<int, int> locateVertex(double x, double y) const;

  static bool isNoTopology(const std::pair<int, int>& topo)
  {
    return topo.first == -1 && topo.second == -1;
  }

  bool isValidTopology(const std::pair<int, int>& topo) const
  {
    if (topo.first < 0 || topo.second < 0)
      return false;
    const size_t obs_index = static_cast<size_t>(topo.first);
    if (obs_index >= obs_.size())
      return false;
    return static_cast<size_t>(topo.second) < obs_[obs_index].ordered_vertices_.size();
  }

  inline bool areConsecutive(const std::pair<int, int>& prev, const std::pair<int, int>& next) const
  {
    if (!isValidTopology(prev) || !isValidTopology(next) || prev.first != next.first)
      return false;
    const auto& vertices = obs_[static_cast<size_t>(prev.first)].ordered_vertices_;
    const size_t expected_next = (static_cast<size_t>(prev.second) + 1) % vertices.size();
    return expected_next == static_cast<size_t>(next.second);
  }

  // Returns the (x, y) coordinate of the obstacle vertex immediately before
  // `curr` in its polygon's `ordered_vertices_` list, with wrap-around
  // (so the predecessor of vertex 0 is the last vertex of that polygon).
  //
  // `curr` is a topological index:
  //   curr.first  -> index into `obs_`  (which obstacle)
  //   curr.second -> index into that obstacle's `ordered_vertices_`
  //
  // NOTE: the return value is the vertex *coordinate* (x, y), NOT a topological
  // index. If topology identity is needed, compute the predecessor/successor
  // index within curr.first directly; a coordinate round-trip through
  // locateVertex() is ambiguous when obstacle contours share a coordinate.
  //
  // Invalid topology is reported with std::nullopt so callers can propagate a
  // normal planning failure instead of allowing an exception to escape.
  inline std::optional<std::pair<int, int>> getPrevObs(
    const std::pair<int, int>& curr) const
  {
    if (!isValidTopology(curr))
      return std::nullopt;
    const auto& vertices = obs_[static_cast<size_t>(curr.first)].ordered_vertices_;
    const size_t index = (static_cast<size_t>(curr.second) + vertices.size() - 1) %
      vertices.size();
    return vertices[index];
  }

  inline std::optional<std::pair<int, int>> getNextObs(
    const std::pair<int, int>& curr) const
  {
    if (!isValidTopology(curr))
      return std::nullopt;
    const auto& vertices = obs_[static_cast<size_t>(curr.first)].ordered_vertices_;
    const size_t index = (static_cast<size_t>(curr.second) + 1) % vertices.size();
    return vertices[index];
  }

  [[nodiscard]] bool calculateVisibilityRegion(int x, int y,
    VisibilityRegion& visibility_region);

  [[nodiscard]] OperationStatus calculateVisibilityRegion(
    int x, int y, VisibilityRegion& visibility_region,
    const StopToken& stop_token);

  // Legacy projection with the same exact-vertex-only topology semantics as
  // the compatibility getVisibilityRegion() overload.
  [[nodiscard]] bool calculateVisibilityRegion(int x, int y,
    std::vector<std::pair<double, double>>& result_V,
    std::vector<std::pair<int, int>>& topo_V);

  bool validateBoundaryEndpoint(const BoundaryEndpoint& endpoint,
    std::string* error = nullptr) const;
  [[nodiscard]] OperationStatus validateBoundaryEndpoint(
    const BoundaryEndpoint& endpoint, const StopToken& stop_token,
    std::string* error = nullptr) const;
  // Source-independent compatibility check for endpoint metadata. Generated
  // full/scoped regions additionally pass the internal source-aware exact
  // geometry validator before they are cached or consumed by the planner.
  bool validateVisibilityRegion(const VisibilityRegion& visibility_region,
    std::string* error = nullptr) const;
  [[nodiscard]] OperationStatus validateVisibilityRegion(
    const VisibilityRegion& visibility_region, const StopToken& stop_token,
    std::string* error = nullptr) const;
  bool boundarySupportsConsecutive(const BoundaryEndpoint& prev,
    const BoundaryEndpoint& next) const;
  [[nodiscard]] OperationStatus boundarySupportsConsecutive(
    const BoundaryEndpoint& prev, const BoundaryEndpoint& next,
    bool& supports, const StopToken& stop_token) const;

  struct CDTEdge
  {
    std::pair<int, int> a, b;
    bool is_constrained;
  };
  // Return at most max_edges so visualization callers cannot first copy an
  // unbounded CDT edge list and only then discover their output budget.
  std::vector<CDTEdge> getCDTEdges(
    size_t max_edges = std::numeric_limits<size_t>::max()) const;

#ifdef RAYSTAR_TESTING
public:  // Legacy white-box tests; production builds keep internals private.
#else
private:
#endif
  // Function-pointer injection is private and exists only so failure paths of
  // CGAL validation can be covered deterministically by unit tests.
  using CdtValidator = bool (*)(
    const constrained_delaunay_triangulation::CDT&);

  friend class PolymapTestPeer;

  Polymap(const GridMap& grid_map, int start_x, int start_y,
    int goal_x, int goal_y, const StopToken& stop_token = StopToken{});
  Polymap(const GridMap& grid_map, int start_x, int start_y,
    int goal_x, int goal_y, const Point2d& start_position,
    const Point2d& goal_position,
    const StopToken& stop_token = StopToken{});

  [[nodiscard]] static PolymapCreateResult finishCreation(
    Polymap&& candidate);

  bool getPolyObstacles(int start_x, int start_y, int goal_x, int goal_y);
  [[nodiscard]] OperationStatus getPolyObstacles(
    int start_x, int start_y, int goal_x, int goal_y,
    const StopToken& stop_token);

  void simplifyPolyObstacles(int start_x, int start_y, int goal_x, int goal_y);
  void simplifyPolyObstacles(const Point2d& start, const Point2d& goal);

  struct VisibilityCacheEntry
  {
    VisibilityRegion vertices;
  };

  std::vector<int> vertices_location_x_flat_;
  std::vector<int> vertices_location_y_flat_;

  std::unordered_map<int, VisibilityCacheEntry> visibility_storage_;

  inline int vertexIndex(int x, int y) const { return y * xsize_ + x; }

  bool calculateVisibilityRegionImpl(int x, int y,
    VisibilityRegion& visibility_region, const StopToken& stop_token);
  bool calculateVisibilityRegionImpl(
    const Point2d& source, const exact_geometry::Point& exact_source,
    const std::optional<std::pair<int, int>>& obstacle_vertex_source,
    VisibilityRegion& visibility_region, const StopToken& stop_token);
  bool getVisibilityRegionImpl(int start_x, int start_y,
    VisibilityRegion& visibility_region, const StopToken& stop_token,
    std::string* error);
  bool validateBoundaryEndpointImpl(const BoundaryEndpoint& endpoint,
    const StopToken& stop_token, std::string* error) const;
  bool validateVisibilityRegionImpl(const VisibilityRegion& visibility_region,
    const StopToken& stop_token, std::string* error) const;
  bool boundarySupportsConsecutiveImpl(const BoundaryEndpoint& prev,
    const BoundaryEndpoint& next, bool& supports,
    const StopToken& stop_token) const;
  bool getPolyObstaclesImpl(int start_x, int start_y, int goal_x, int goal_y,
    const StopToken& stop_token);
  bool constructCGALRelatedImpl(
    CdtValidator validator, std::string& error,
    const StopToken& stop_token);
  bool validateObstacleTopologyImpl(
    std::string& error, bool validate_edge_relations,
    const StopToken& stop_token) const;
  bool simplificationChordIsTopologicallySafeImpl(
    size_t obstacle_index, size_t current_index,
    const StopToken& stop_token) const;
  bool registerVerticesImpl(
    std::string& error, const StopToken& stop_token);
  bool simplifyPolyObstaclesImpl(
    const Point2d& start, const Point2d& goal,
    const StopToken& stop_token);
  bool validateFreeSpaceInteriorImpl(
    const Point2d& point, const StopToken& stop_token,
    std::string* error) const;

  bool isInTri(int x1, int y1, int x2, int y2, int x3, int y3, double x, double y);

  BoundaryEndpoint makeEndpoint(const Point2d& position) const;
  BoundaryEndpoint makeEndpoint(const std::pair<int, int>& position) const;
  BoundaryEndpoint makeIntersectionEndpoint(const exact_geometry::Point& position,
    const DirectedObstacleEdge& supporting_edge) const;
  std::optional<DirectedObstacleEdge> obstacleEdgeBetween(
    const Point2d& from, const Point2d& to) const;
  bool isAnObstacleEdge(const Point2d& from, const Point2d& to) const;
  bool pointLiesOnObstacleEdge(const BoundaryEndpoint& point,
    const DirectedObstacleEdge& edge) const;

  static std::pair<int, int> toLegacyTopology(const ObstacleVertexId& vertex)
  {
    return {vertex.obstacle, vertex.vertex};
  }

  static ObstacleVertexId fromLegacyTopology(const std::pair<int, int>& vertex)
  {
    return {vertex.first, vertex.second};
  }

  [[nodiscard]] bool constructCGALRelated(
    CdtValidator validator, std::string& error);
  [[nodiscard]] OperationStatus constructCGALRelated(
    CdtValidator validator, std::string& error,
    const StopToken& stop_token);
  void clearCGALRelatedState();
  void clearStoppedConstructionState();
  [[nodiscard]] bool validateObstacleTopology(
    std::string& error, bool validate_edge_relations) const;
  [[nodiscard]] OperationStatus validateObstacleTopology(
    std::string& error, bool validate_edge_relations,
    const StopToken& stop_token) const;
  [[nodiscard]] bool simplificationChordIsTopologicallySafe(
    size_t obstacle_index, size_t current_index) const;
  [[nodiscard]] OperationStatus simplificationChordIsTopologicallySafe(
    size_t obstacle_index, size_t current_index, bool& safe,
    const StopToken& stop_token) const;
  [[nodiscard]] bool registerVertices(std::string& error);
  [[nodiscard]] OperationStatus registerVertices(
    std::string& error, const StopToken& stop_token);
  [[nodiscard]] OperationStatus simplifyPolyObstacles(
    int start_x, int start_y, int goal_x, int goal_y,
    const StopToken& stop_token);
  [[nodiscard]] OperationStatus simplifyPolyObstacles(
    const Point2d& start, const Point2d& goal,
    const StopToken& stop_token);

  bool isFacetInsideObstacle(int facet_idx) const;

  int xsize_ = 0;
  int ysize_ = 0;
  std::vector<uint8_t> data_;
  std::vector<Obs> obs_;

  // Grid flood-fill connectivity only; CDT construction has an independent
  // status so a geometry failure is not misreported as an unreachable goal.
  bool solution_exist_ = false;
  bool no_path_ = false;

  bool cdt_ready_ = false;
  bool construction_stopped_ = false;
  std::string construction_error_;
  constrained_delaunay_triangulation::CDT cdt_;
  std::unordered_map<long long, int> cdt_table_;
  int cdt_ver_num_ = 0;
  std::vector<std::vector<std::pair<int, int>>> facets_;

  inline int locateAdjacentFacet(std::pair<int, int> prev, std::pair<int, int> next) const
  {
    auto it = cdt_table_.find(static_cast<long long>(prev.first + prev.second * xsize_) +
      static_cast<long long>(next.first + next.second * xsize_) * static_cast<long long>(xsize_ * ysize_));
    if (it == cdt_table_.end()) return -1;
    return it->second;
  }
};

struct PolymapCreateResult
{
  PolymapCreateStatus status = PolymapCreateStatus::failure;
  std::optional<Polymap> value;
  std::string error;

  [[nodiscard]] explicit operator bool() const noexcept
  {
    return status == PolymapCreateStatus::ready && value.has_value();
  }
};

}  // namespace raystar
