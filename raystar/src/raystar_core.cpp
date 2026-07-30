#include <raystar/raystar_core.h>
#include <chrono>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <limits>
#include <exception>
#include <optional>

#include "exact_point_location.h"
#include "visibility_validation.h"

namespace raystar {

bool validateMapResourceBudget(size_t width,
                               size_t height,
                               size_t data_size,
                               const PlanningLimits& limits,
                               MapResourceEstimate& estimate,
                               std::string& error) {
  estimate = MapResourceEstimate{};
  if (width == 0 || height == 0) {
    error = "Invalid map: width and height must be positive";
    return false;
  }
  if (width > std::numeric_limits<size_t>::max() / height) {
    error = "Invalid map: width * height overflows size_t";
    return false;
  }

  const size_t cell_count = width * height;
  const size_t max_int = static_cast<size_t>(std::numeric_limits<int>::max());
  if (width > max_int || height > max_int) {
    error = "Invalid map: width and height must fit signed int indexing";
    return false;
  }
  if (cell_count > max_int) {
    error = "Invalid map: cell count must fit signed int indexing";
    return false;
  }
  if (data_size != cell_count) {
    error = "Invalid map: data size is " + std::to_string(data_size) + ", expected " +
            std::to_string(cell_count);
    return false;
  }
  if (limits.max_map_cells == 0) {
    error = "Invalid planning limits: max_map_cells must be greater than zero";
    return false;
  }
  if (limits.max_map_bytes == 0) {
    error = "Invalid planning limits: max_map_bytes must be greater than zero";
    return false;
  }
  if (limits.max_map_bytes < kEstimatedPlannerMapBytesPerCell) {
    error = "Invalid planning limits: max_map_bytes must be at least " +
            std::to_string(kEstimatedPlannerMapBytesPerCell);
    return false;
  }
  if (cell_count > limits.max_map_cells) {
    error = "Invalid map: cell count " + std::to_string(cell_count) +
            " exceeds max_map_cells=" + std::to_string(limits.max_map_cells);
    return false;
  }
  if (cell_count > std::numeric_limits<size_t>::max() / kEstimatedPlannerMapBytesPerCell) {
    error = "Invalid map: estimated planner map memory overflows size_t";
    return false;
  }

  const size_t estimated_bytes = cell_count * kEstimatedPlannerMapBytesPerCell;
  if (estimated_bytes > limits.max_map_bytes) {
    error = "Invalid map: estimated planner map memory " + std::to_string(estimated_bytes) +
            " bytes exceeds max_map_bytes=" + std::to_string(limits.max_map_bytes);
    return false;
  }

  estimate.cell_count = cell_count;
  estimate.estimated_bytes = estimated_bytes;
  error.clear();
  return true;
}

namespace {

bool isRepresentableAsInt(double value) {
  return std::isfinite(value) && value >= static_cast<double>(std::numeric_limits<int>::lowest()) &&
         value <= static_cast<double>(std::numeric_limits<int>::max());
}

int safeIntCast(double value) {
  return isRepresentableAsInt(value) ? static_cast<int>(value) : 0;
}

int safeRoundedIntCast(double value) {
  if (!isRepresentableAsInt(value))
    return 0;
  const double rounded = std::round(value);
  return isRepresentableAsInt(rounded) ? static_cast<int>(rounded) : 0;
}

// Validate the map shape before any caller-owned GridMap is copied into the
// planner's working map.  Besides making the multiplication overflow-safe,
// this is the single admission point for the configurable cell budget.  The
// signed-indexing check remains an independent invariant because several
// legacy geometry tables use int coordinates even when the configured budget
// is larger.
bool validateMapShapeAndAdmission(const GridMap& grid_map,
                                  const PlanningLimits& limits,
                                  size_t& width,
                                  size_t& height,
                                  size_t& cell_count,
                                  std::string& error) {
  width = static_cast<size_t>(grid_map.width);
  height = static_cast<size_t>(grid_map.height);
  MapResourceEstimate estimate;
  if (!validateMapResourceBudget(width, height, grid_map.data.size(), limits, estimate, error)) {
    return false;
  }
  cell_count = estimate.cell_count;
  return true;
}

bool validateNonMapPlanningLimits(const PlanningLimits& limits, std::string& error) {
  const size_t max_int = static_cast<size_t>(std::numeric_limits<int>::max());
  if (limits.max_k <= 0) {
    error = "Invalid planning limits: max_k must be greater than zero";
    return false;
  }
  if (limits.max_nodes == 0 || limits.max_nodes > max_int) {
    error = "Invalid planning limits: max_nodes must be between 1 and " + std::to_string(max_int);
    return false;
  }
  if (limits.planning_timeout.count() < 0) {
    error = "Invalid planning limits: planning timeout must not be negative";
    return false;
  }
  if (limits.max_path_points < 2) {
    error = "Invalid planning limits: max_path_points must be at least 2";
    return false;
  }
  // Zero is an intentional opt-out for the ROS debug arrays/markers.  It
  // never affects the search-node admission limit itself.
  if (limits.max_debug_nodes > max_int) {
    error =
      "Invalid planning limits: max_debug_nodes must be between 0 and " + std::to_string(max_int);
    return false;
  }
  if (limits.max_response_bytes == 0) {
    error = "Invalid planning limits: max_response_bytes must be greater than zero";
    return false;
  }
  return true;
}

std::string boundedPlanningExceptionMessage(const char* detail) {
  constexpr size_t kMaxMessageBytes = 1024;
  constexpr const char* prefix = "Planning failed with exception: ";
  std::string message(prefix);
  if (message.size() >= kMaxMessageBytes || detail == nullptr) {
    message.resize(std::min(message.size(), kMaxMessageBytes));
    return message;
  }
  const size_t remaining = kMaxMessageBytes - message.size();
  size_t detail_bytes = 0;
  while (detail_bytes < remaining && detail[detail_bytes] != '\0') ++detail_bytes;
  message.append(detail, detail_bytes);
  return message;
}

bool validatePlanRequest(const GridMap& grid_map,
                         int start_x,
                         int start_y,
                         int goal_x,
                         int goal_y,
                         int K,
                         const PlanningLimits& limits,
                         std::string& error) {
  size_t width = 0;
  size_t height = 0;
  size_t cell_count = 0;
  if (!validateMapShapeAndAdmission(grid_map, limits, width, height, cell_count, error)) {
    return false;
  }
  (void)cell_count;
  if (!validateNonMapPlanningLimits(limits, error))
    return false;

  if (K <= 0) {
    error = "Invalid K: must be greater than zero";
    return false;
  }
  if (K > limits.max_k) {
    error = "Invalid K: requested " + std::to_string(K) +
            " exceeds max_k=" + std::to_string(limits.max_k);
    return false;
  }

  const auto is_in_bounds = [width, height](int x, int y) {
    return x >= 0 && y >= 0 && static_cast<size_t>(x) < width && static_cast<size_t>(y) < height;
  };

  if (!is_in_bounds(start_x, start_y)) {
    error = "Invalid start: coordinates (" + std::to_string(start_x) + ", " +
            std::to_string(start_y) + ") are outside map bounds";
    return false;
  }
  if (!is_in_bounds(goal_x, goal_y)) {
    error = "Invalid goal: coordinates (" + std::to_string(goal_x) + ", " + std::to_string(goal_y) +
            ") are outside map bounds";
    return false;
  }

  error.clear();
  return true;
}

std::string planningLimitMessage(PlanningLimitReached limit,
                                 const PlanningLimits& limits,
                                 size_t path_count) {
  std::string message;
  if (limit == PlanningLimitReached::max_nodes) {
    message = "Planning stopped after reaching max_nodes=" + std::to_string(limits.max_nodes);
  } else if (limit == PlanningLimitReached::max_path_points) {
    message =
      "Planning stopped before exceeding max_path_points=" + std::to_string(limits.max_path_points);
  } else if (limit == PlanningLimitReached::timeout) {
    message = "Planning stopped after reaching planning timeout of " +
              std::to_string(limits.planning_timeout.count()) + " ms";
  } else if (limit == PlanningLimitReached::cancelled) {
    message = "Planning was cancelled";
  } else {
    message = "Planning stopped";
  }

  if (path_count == 0) {
    message += limit == PlanningLimitReached::cancelled ? "; no path was found before cancellation"
                                                        : "; no path was found before the limit";
  } else {
    message +=
      "; returning " + std::to_string(path_count) + " path(s), but the search is incomplete";
  }
  return message;
}

std::pair<int, int> legacyTopology(const BoundaryEndpoint& endpoint) {
  if (const auto* vertex = std::get_if<ObstacleVertexId>(&endpoint.support))
    return {vertex->obstacle, vertex->vertex};
  return {-1, -1};
}

bool projectVisibilityRegion(const VisibilityRegion& visibility_region,
                             std::vector<std::pair<double, double>>& positions,
                             std::vector<std::pair<int, int>>& topology,
                             const StopToken& stop_token) {
  positions.clear();
  topology.clear();
  if (stop_token.poll())
    return false;
  positions.reserve(visibility_region.size());
  topology.reserve(visibility_region.size());
  for (const auto& endpoint : visibility_region) {
    if (stop_token.poll()) {
      positions.clear();
      topology.clear();
      return false;
    }
    positions.emplace_back(endpoint.position);
    topology.emplace_back(legacyTopology(endpoint));
  }
  if (stop_token.poll()) {
    positions.clear();
    topology.clear();
    return false;
  }
  return true;
}

BoundaryEndpoint exactEndpoint(const std::pair<int, int>& position,
                               const std::pair<int, int>& topology) {
  BoundaryEndpoint endpoint;
  endpoint.position = {static_cast<double>(position.first), static_cast<double>(position.second)};
  endpoint.exact_position = exact_geometry::Point(position.first, position.second);
  endpoint.support = ObstacleVertexId{topology.first, topology.second};
  return endpoint;
}

bool canonicalizeExactObstacleVertices(const Polymap* pMap,
                                       VisibilityRegion& visibility_region,
                                       const StopToken& stop_token) {
  if (!pMap)
    return true;
  if (stop_token.poll())
    return false;
  constexpr double coordinate_tolerance = 1e-9;
  for (auto& endpoint : visibility_region) {
    if (stop_token.poll())
      return false;
    const auto* vertex_id = std::get_if<ObstacleVertexId>(&endpoint.support);
    if (!vertex_id)
      continue;
    const std::pair<int, int> topology = {vertex_id->obstacle, vertex_id->vertex};
    if (!pMap->isValidTopology(topology))
      continue;
    const auto& vertex = pMap->obstacles()[static_cast<size_t>(vertex_id->obstacle)]
                           .ordered_vertices_[static_cast<size_t>(vertex_id->vertex)];
    if (!std::isfinite(endpoint.position.first) || !std::isfinite(endpoint.position.second) ||
        std::abs(endpoint.position.first - static_cast<double>(vertex.first)) >
          coordinate_tolerance ||
        std::abs(endpoint.position.second - static_cast<double>(vertex.second)) >
          coordinate_tolerance) {
      continue;
    }
    endpoint.position = {static_cast<double>(vertex.first), static_cast<double>(vertex.second)};
    endpoint.exact_position = exact_geometry::Point(vertex.first, vertex.second);
  }
  return true;
}

VisibilityRegion adaptLegacyVisibilityRegion(
  const Polymap* pMap,
  const std::vector<std::pair<double, double>>& positions,
  const std::vector<std::pair<int, int>>& topology,
  bool& valid,
  std::string& error) {
  VisibilityRegion result;
  valid = true;
  error.clear();

  if (!pMap) {
    valid = false;
    error = "Cannot adapt a legacy visibility region without a Polymap";
    return result;
  }
  if (positions.size() != topology.size()) {
    valid = false;
    error = "Visibility region and topology sizes do not match";
    return result;
  }

  result.reserve(positions.size());
  for (size_t i = 0; i < positions.size(); ++i) {
    BoundaryEndpoint endpoint;
    endpoint.position = positions[i];
    if (!isRepresentableAsInt(endpoint.position.first) ||
        !isRepresentableAsInt(endpoint.position.second)) {
      valid = false;
      error =
        "Visibility region contains a vertex outside the representable int range "
        "at position " +
        std::to_string(i);
      return {};
    }
    endpoint.exact_position =
      exact_geometry::Point(endpoint.position.first, endpoint.position.second);

    const auto& legacy = topology[i];
    if (Polymap::isNoTopology(legacy)) {
      endpoint.support = std::monostate{};
    } else if (!pMap->isValidTopology(legacy)) {
      endpoint.support = ObstacleVertexId{legacy.first, legacy.second};
    } else {
      const auto& obstacle = pMap->obstacles()[static_cast<size_t>(legacy.first)];
      const auto& vertex = obstacle.ordered_vertices_[static_cast<size_t>(legacy.second)];
      constexpr double coordinate_tolerance = 1e-9;
      if (std::abs(endpoint.position.first - static_cast<double>(vertex.first)) <=
            coordinate_tolerance &&
          std::abs(endpoint.position.second - static_cast<double>(vertex.second)) <=
            coordinate_tolerance) {
        endpoint.position = {static_cast<double>(vertex.first), static_cast<double>(vertex.second)};
        endpoint.exact_position = exact_geometry::Point(vertex.first, vertex.second);
        endpoint.support = ObstacleVertexId{legacy.first, legacy.second};
      } else {
        const int next_vertex =
          (legacy.second + 1) % static_cast<int>(obstacle.ordered_vertices_.size());
        endpoint.support = DirectedObstacleEdge{ObstacleVertexId{legacy.first, legacy.second},
                                                ObstacleVertexId{legacy.first, next_vertex}};
      }
    }
    result.emplace_back(std::move(endpoint));
  }
  return result;
}

bool validatePlanEndpointShape(const GridMap& grid_map,
                               const PlanEndpoint& endpoint,
                               const char* label,
                               std::string& error) {
  // Check the geometric value before deriving diagnostics from `cell_`.
  // The Point2d overload intentionally uses {-1,-1} as a sentinel when it
  // receives NaN/Inf; reporting an out-of-bounds cell in that case hides the
  // actual input error.
  if (!std::isfinite(endpoint.position_.first) || !std::isfinite(endpoint.position_.second)) {
    error = std::string("Invalid ") + label + ": coordinates must be finite";
    return false;
  }

  const auto& cell = endpoint.cell_;
  if (cell.first < 0 || cell.second < 0 || static_cast<size_t>(cell.first) >= grid_map.width ||
      static_cast<size_t>(cell.second) >= grid_map.height) {
    error = std::string("Invalid ") + label + ": cell coordinates (" + std::to_string(cell.first) +
            ", " + std::to_string(cell.second) + ") are outside map bounds";
    return false;
  }
  const double floor_x = std::floor(endpoint.position_.first);
  const double floor_y = std::floor(endpoint.position_.second);
  if (floor_x != static_cast<double>(cell.first) || floor_y != static_cast<double>(cell.second)) {
    error = std::string("Invalid ") + label + ": continuous position and grid cell disagree";
    return false;
  }
  const auto index =
    static_cast<size_t>(cell.second) * grid_map.width + static_cast<size_t>(cell.first);
  if (index >= grid_map.data.size() || grid_map.data[index] != 0) {
    error = std::string("Invalid ") + label + ": corresponding occupancy cell is occupied";
    return false;
  }
  return true;
}

bool validatePlanRequest(const GridMap& grid_map,
                         const PlanEndpoint& start,
                         const PlanEndpoint& goal,
                         int K,
                         const PlanningLimits& limits,
                         std::string& error) {
  // Keep all map/limit/K diagnostics byte-for-byte compatible with the legacy
  // request validator, then apply the endpoint-specific contract.
  size_t width = 0;
  size_t height = 0;
  size_t cell_count = 0;
  if (!validateMapShapeAndAdmission(grid_map, limits, width, height, cell_count, error)) {
    return false;
  }
  (void)cell_count;
  if (!validateNonMapPlanningLimits(limits, error))
    return false;
  if (K <= 0) {
    error = "Invalid K: must be greater than zero";
    return false;
  }
  if (K > limits.max_k) {
    error = "Invalid K: requested " + std::to_string(K) +
            " exceeds max_k=" + std::to_string(limits.max_k);
    return false;
  }
  if (!validatePlanEndpointShape(grid_map, start, "start", error) ||
      !validatePlanEndpointShape(grid_map, goal, "goal", error)) {
    return false;
  }
  error.clear();
  return true;
}

OperationStatus hybridPathSelfCrosses(const Point2d& start,
                                      const std::vector<std::pair<int, int>>& turns,
                                      const Point2d& candidate,
                                      bool& self_crosses,
                                      const StopToken& stop_token) {
  self_crosses = false;
  if (stop_token.poll())
    return OperationStatus::stopped;
  std::vector<exact_geometry::Point> points;
  points.reserve(turns.size() + 2);
  points.emplace_back(start.first, start.second);
  for (const auto& turn : turns) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    points.emplace_back(turn.first, turn.second);
  }
  points.emplace_back(candidate.first, candidate.second);

  const auto& new_point = points.back();
  for (size_t i = 0; i + 1 < points.size() - 1; ++i) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    if (points[i] == new_point) {
      self_crosses = true;
      return OperationStatus::success;
    }
  }
  if (points.size() >= 4) {
    const auto new_segment =
      exact_geometry::Kernel::Segment_2(points[points.size() - 2], points.back());
    for (size_t i = 0; i + 2 < points.size() - 1; ++i) {
      if (stop_token.poll())
        return OperationStatus::stopped;
      const auto old_segment = exact_geometry::Kernel::Segment_2(points[i], points[i + 1]);
      if (CGAL::do_intersect(old_segment, new_segment)) {
        self_crosses = true;
        return OperationStatus::success;
      }
    }
  }
  return OperationStatus::success;
}

OperationStatus generateChildrenFromSource(const Polymap* pMap,
                                           int owner_index,
                                           const Point2d& source,
                                           const exact_geometry::Point& exact_source,
                                           double start_angle,
                                           double end_angle,
                                           double gcost,
                                           const VisibilityRegion& visibility_region,
                                           std::vector<Child>& output,
                                           const StopToken& stop_token,
                                           std::string* error) {
  output.clear();
  if (error)
    error->clear();
  const auto fail = [&](const std::string& message) {
    output.clear();
    if (error)
      *error = message;
    return OperationStatus::failure;
  };
  const auto stopped = [&]() {
    output.clear();
    if (error)
      error->clear();
    return OperationStatus::stopped;
  };
  if (stop_token.poll())
    return stopped();
  if (!pMap)
    return fail("Cannot generate children without a Polymap");
  if (visibility_region.size() < 2)
    return fail("Visibility region has fewer than two vertices");

  std::string validation_error;
  const auto visibility_status =
    pMap->validateVisibilityRegion(visibility_region, stop_token, &validation_error);
  if (visibility_status == OperationStatus::stopped)
    return stopped();
  if (visibility_status == OperationStatus::failure)
    return fail(validation_error);

  std::vector<double> theta_list(visibility_region.size());
  for (size_t i = 0; i < visibility_region.size(); ++i) {
    if (stop_token.poll())
      return stopped();
    const auto& position = visibility_region[i].position;
    theta_list[i] = atan2(position.second - source.second, position.first - source.first);
    theta_list[i] = start_angle + normalize_angle_positive(theta_list[i] - start_angle);
  }

  std::vector<bool> is_a_left_gap(visibility_region.size(), false);
  std::vector<size_t> valid_child_indices;
  for (size_t i = 0; i + 1 < visibility_region.size(); ++i) {
    if (stop_token.poll())
      return stopped();
    const size_t next = i + 1;
    const auto& current_endpoint = visibility_region[i];
    const auto& next_endpoint = visibility_region[next];
    if (exactPoint(current_endpoint) == exactPoint(next_endpoint))
      continue;
    if (exact_geometry::isSameDirectedRay(
          exact_source, exactPoint(current_endpoint), exactPoint(next_endpoint))) {
      if (std::holds_alternative<std::monostate>(current_endpoint.support) &&
          std::holds_alternative<std::monostate>(next_endpoint.support)) {
        continue;
      }
      bool supports_consecutive = false;
      const auto support_status = pMap->boundarySupportsConsecutive(
        next_endpoint, current_endpoint, supports_consecutive, stop_token);
      if (support_status == OperationStatus::stopped)
        return stopped();
      if (support_status == OperationStatus::failure)
        return fail("Visibility boundary support validation failed");
      if (!supports_consecutive) {
        const auto current_distance =
          (exactPoint(current_endpoint) - exact_source).squared_length();
        const auto next_distance = (exactPoint(next_endpoint) - exact_source).squared_length();
        if (current_distance > next_distance)
          is_a_left_gap[i] = true;
        valid_child_indices.emplace_back(i);
      }
    }
  }

  struct PendingChild {
    BoundaryEndpoint near_endpoint;
    BoundaryEndpoint far_endpoint;
    bool is_left_gap = false;
    double start_angle = 0.0;
    double end_angle = 0.0;
  };
  std::vector<PendingChild> pending_children;
  pending_children.reserve(valid_child_indices.size());
  for (const auto i : valid_child_indices) {
    if (stop_token.poll())
      return stopped();
    const size_t next = i + 1;
    const size_t near_index = is_a_left_gap[i] ? next : i;
    const size_t far_index = is_a_left_gap[i] ? i : next;
    const auto& near_endpoint = visibility_region[near_index];
    const auto& far_endpoint = visibility_region[far_index];
    const auto* near_vertex = std::get_if<ObstacleVertexId>(&near_endpoint.support);
    if (!near_vertex) {
      return fail("Visibility gap has an edge-interior or unsupported near anchor at position " +
                  std::to_string(near_index));
    }
    const std::pair<int, int> near_topology = {near_vertex->obstacle, near_vertex->vertex};
    if (!pMap->isValidTopology(near_topology))
      return fail("Visibility gap has an invalid near anchor topology");

    double child_start_angle = 0.0;
    double child_end_angle = 0.0;
    if (is_a_left_gap[i]) {
      const auto next_obs = pMap->getNextObs(near_topology);
      if (!next_obs)
        return fail("Cannot find the next obstacle vertex for a visibility gap");
      child_start_angle = normalize_angle(theta_list[next]);
      const double contour_angle = atan2(next_obs->second - near_endpoint.position.second,
                                         next_obs->first - near_endpoint.position.first);
      child_end_angle =
        child_start_angle + normalize_angle_positive(contour_angle - child_start_angle);
    } else {
      const auto prev_obs = pMap->getPrevObs(near_topology);
      if (!prev_obs)
        return fail("Cannot find the previous obstacle vertex for a visibility gap");
      const double contour_angle = atan2(prev_obs->second - near_endpoint.position.second,
                                         prev_obs->first - near_endpoint.position.first);
      child_start_angle = contour_angle;
      child_end_angle = contour_angle + normalize_angle_positive(theta_list[i] - contour_angle);
    }

    // The near endpoint is topologically an obstacle vertex.  Recover the
    // integer coordinate from that identity instead of rounding a projected
    // double, so the Child contract remains strictly integral even for a
    // continuous root source.
    const auto& canonical = pMap->obstacles()[static_cast<size_t>(near_topology.first)]
                              .ordered_vertices_[static_cast<size_t>(near_topology.second)];
    BoundaryEndpoint canonical_near_endpoint = near_endpoint;
    canonical_near_endpoint.position = {static_cast<double>(canonical.first),
                                        static_cast<double>(canonical.second)};
    canonical_near_endpoint.exact_position =
      exact_geometry::Point(canonical.first, canonical.second);
    pending_children.push_back(PendingChild{std::move(canonical_near_endpoint),
                                            far_endpoint,
                                            is_a_left_gap[i],
                                            child_start_angle,
                                            child_end_angle});
  }

  std::vector<Child> generated_children;
  generated_children.reserve(pending_children.size());
  for (size_t i = 0; i < pending_children.size(); ++i) {
    if (stop_token.poll())
      return stopped();
    auto& pending = pending_children[i];
    const double child_gcost =
      gcost + std::hypot(source.first - pending.near_endpoint.position.first,
                         source.second - pending.near_endpoint.position.second);
    generated_children.emplace_back(owner_index,
                                    static_cast<int>(i),
                                    pending.near_endpoint,
                                    pending.far_endpoint,
                                    pending.is_left_gap,
                                    pending.start_angle,
                                    pending.end_angle,
                                    child_gcost);
  }
  if (stop_token.poll())
    return stopped();
  output = std::move(generated_children);
  (void)end_angle;  // retained in the helper signature for the Node contract
  return OperationStatus::success;
}

}  // namespace

static detail::VisibilityBoundaryMode visibilityBoundaryModeForNode(const Polymap& polymap,
                                                                    const Node& node) {
  if (node.isContinuousRoot())
    return detail::VisibilityBoundaryMode::closed_cycle;
  if (node.parentIndex() >= 0)
    return detail::VisibilityBoundaryMode::open_sector;
  if (node.hasValidSeed() &&
      polymap.isValidTopology(polymap.locateVertex(node.seed().first, node.seed().second))) {
    return detail::VisibilityBoundaryMode::open_sector;
  }
  return detail::VisibilityBoundaryMode::closed_cycle;
}

Child::Child(int nindex,
             int cindex,
             const BoundaryEndpoint& c_endpoint,
             const BoundaryEndpoint& o_endpoint,
             bool is_left_gap)
  : Nindex_(nindex)
  , Cindex_(cindex)
  , start_angle_(0)
  , end_angle_(0)
  , c_({safeRoundedIntCast(c_endpoint.position.first),
        safeRoundedIntCast(c_endpoint.position.second)})
  , o_(o_endpoint.position)
  , c_endpoint_(c_endpoint)
  , o_endpoint_(o_endpoint)
  , c_obs_index_(-1)
  , c_ver_index_(-1)
  , o_obs_index_(-1)
  , o_ver_index_(-1)
  , is_a_left_gap_(is_left_gap)
  , c_gcost_(0)
  , c_hcost_(0) {
  const auto c_topology = legacyTopology(c_endpoint_);
  c_obs_index_ = c_topology.first;
  c_ver_index_ = c_topology.second;
  const auto o_topology = legacyTopology(o_endpoint_);
  o_obs_index_ = o_topology.first;
  o_ver_index_ = o_topology.second;
}

Child::Child(int nindex,
             int cindex,
             const BoundaryEndpoint& c_endpoint,
             const BoundaryEndpoint& o_endpoint,
             bool is_left_gap,
             double start_angle,
             double end_angle,
             double g_cost)
  : Child(nindex, cindex, c_endpoint, o_endpoint, is_left_gap) {
  start_angle_ = start_angle;
  end_angle_ = end_angle;
  c_gcost_ = g_cost;
}

Node::Node(const Polymap* pMap,
           int Nindex,
           double start_x,
           double start_y,
           double Gcost,
           double Hcost,
           const VisibilityRegion& visibility_region)
  : Node(pMap, Nindex, start_x, start_y, Gcost, Hcost, visibility_region, StopToken{}) {}

Node::Node(const Polymap* pMap,
           int Nindex,
           double start_x,
           double start_y,
           double Gcost,
           double Hcost,
           const VisibilityRegion& visibility_region,
           const StopToken& stop_token)
  : Nindex_(Nindex)
  , seed_({safeIntCast(start_x), safeIntCast(start_y)})
  , seed_is_valid_(isRepresentableAsInt(start_x) && isRepresentableAsInt(start_y))
  , start_angle_(0)
  , end_angle_(kTwoPi)
  , parent_index_(-1)
  , start_o_({0.0, 0.0})
  , end_o_({0.0, 0.0})
  , as_a_child_left_gap_(false)
  , Gcost_(Gcost)
  , Hcost_(Hcost)
  , visibility_region_valid_(true) {
  const auto stop_initialization = [&]() {
    visibility_region_valid_ = false;
    visibility_region_error_.clear();
    visibility_region_.clear();
    V_.clear();
    topo_V_.clear();
  };

  local_shortest_path_.emplace_back(seed_);
  path_node_index_.emplace_back(Nindex_);
  if (stop_token.poll()) {
    stop_initialization();
    return;
  }
  visibility_region_ = visibility_region;
  if (stop_token.poll() ||
      !canonicalizeExactObstacleVertices(pMap, visibility_region_, stop_token) ||
      !projectVisibilityRegion(visibility_region_, V_, topo_V_, stop_token)) {
    stop_initialization();
    return;
  }
  if (!pMap) {
    visibility_region_valid_ = false;
    visibility_region_error_ = "Cannot construct a Node without a Polymap";
  } else {
    const auto status =
      pMap->validateVisibilityRegion(visibility_region_, stop_token, &visibility_region_error_);
    if (status == OperationStatus::stopped)
      stop_initialization();
    else if (status == OperationStatus::failure)
      visibility_region_valid_ = false;
  }
}

Node::Node(const Polymap* pMap,
           int Nindex,
           double seed_x,
           double seed_y,
           double Gcost,
           double Hcost,
           int parent_index,
           const VisibilityRegion& visibility_region)
  : Node(pMap, Nindex, seed_x, seed_y, Gcost, Hcost, parent_index, visibility_region, StopToken{}) {
}

Node::Node(const Polymap* pMap,
           int Nindex,
           double seed_x,
           double seed_y,
           double Gcost,
           double Hcost,
           int parent_index,
           const VisibilityRegion& visibility_region,
           const StopToken& stop_token)
  : Nindex_(Nindex)
  , seed_({safeIntCast(seed_x), safeIntCast(seed_y)})
  , seed_is_valid_(isRepresentableAsInt(seed_x) && isRepresentableAsInt(seed_y))
  , start_angle_(0)
  , end_angle_(0)
  , parent_index_(parent_index)
  , start_o_({0.0, 0.0})
  , end_o_({0.0, 0.0})
  , as_a_child_left_gap_(false)
  , Gcost_(Gcost)
  , Hcost_(Hcost)
  , visibility_region_valid_(true) {
  const auto stop_initialization = [&]() {
    visibility_region_valid_ = false;
    visibility_region_error_.clear();
    visibility_region_.clear();
    V_.clear();
    topo_V_.clear();
  };

  if (stop_token.poll()) {
    stop_initialization();
    return;
  }
  visibility_region_ = visibility_region;
  if (stop_token.poll() ||
      !canonicalizeExactObstacleVertices(pMap, visibility_region_, stop_token) ||
      !projectVisibilityRegion(visibility_region_, V_, topo_V_, stop_token)) {
    stop_initialization();
    return;
  }
  if (!pMap) {
    visibility_region_valid_ = false;
    visibility_region_error_ = "Cannot construct a Node without a Polymap";
  } else {
    const auto status =
      pMap->validateVisibilityRegion(visibility_region_, stop_token, &visibility_region_error_);
    if (status == OperationStatus::stopped)
      stop_initialization();
    else if (status == OperationStatus::failure)
      visibility_region_valid_ = false;
  }
}

Node::Node(const Polymap* pMap,
           int Nindex,
           double start_x,
           double start_y,
           double Gcost,
           double Hcost,
           const std::vector<std::pair<double, double>>& visibility_region,
           const std::vector<std::pair<int, int>>& topo_V)
  : Nindex_(Nindex)
  , seed_({safeIntCast(start_x), safeIntCast(start_y)})
  , seed_is_valid_(isRepresentableAsInt(start_x) && isRepresentableAsInt(start_y))
  , start_angle_(0)
  , end_angle_(kTwoPi)
  , parent_index_(-1)
  , start_o_({0.0, 0.0})
  , end_o_({0.0, 0.0})
  , as_a_child_left_gap_(false)
  , Gcost_(Gcost)
  , Hcost_(Hcost)
  , visibility_region_valid_(true) {
  local_shortest_path_.emplace_back(seed_);
  path_node_index_.emplace_back(Nindex_);
  visibility_region_ = adaptLegacyVisibilityRegion(
    pMap, visibility_region, topo_V, visibility_region_valid_, visibility_region_error_);
  const StopToken no_stop;
  (void)projectVisibilityRegion(visibility_region_, V_, topo_V_, no_stop);
}

Node::Node(const Polymap* pMap,
           int Nindex,
           double seed_x,
           double seed_y,
           double Gcost,
           double Hcost,
           int parent_index,
           const std::vector<std::pair<double, double>>& visibility_region,
           const std::vector<std::pair<int, int>>& topo_V)
  : Nindex_(Nindex)
  , seed_({safeIntCast(seed_x), safeIntCast(seed_y)})
  , seed_is_valid_(isRepresentableAsInt(seed_x) && isRepresentableAsInt(seed_y))
  , start_angle_(0)
  , end_angle_(0)
  , parent_index_(parent_index)
  , start_o_({0.0, 0.0})
  , end_o_({0.0, 0.0})
  , as_a_child_left_gap_(false)
  , Gcost_(Gcost)
  , Hcost_(Hcost)
  , visibility_region_valid_(true) {
  visibility_region_ = adaptLegacyVisibilityRegion(
    pMap, visibility_region, topo_V, visibility_region_valid_, visibility_region_error_);
  const StopToken no_stop;
  (void)projectVisibilityRegion(visibility_region_, V_, topo_V_, no_stop);
}

void Node::setFullVisibilityRegion(const VisibilityRegion& visibility_region) {
  const StopToken no_stop;
  (void)setFullVisibilityRegion(visibility_region, no_stop);
}

OperationStatus Node::setFullVisibilityRegion(const VisibilityRegion& visibility_region,
                                              const StopToken& stop_token) {
  if (stop_token.poll())
    return OperationStatus::stopped;

  VisibilityRegion candidate_region = visibility_region;
  if (stop_token.poll())
    return OperationStatus::stopped;

  std::vector<std::pair<double, double>> candidate_positions;
  candidate_positions.reserve(candidate_region.size());
  for (const auto& endpoint : candidate_region) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    candidate_positions.emplace_back(endpoint.position);
  }

  if (stop_token.poll())
    return OperationStatus::stopped;

  full_visibility_region_ = std::move(candidate_region);
  full_V_ = std::move(candidate_positions);
  return OperationStatus::success;
}

bool Node::generateChild(const Polymap* pMap) {
  return generateChild(pMap, nullptr);
}

bool Node::generateChild(const Polymap* pMap, std::string* error) {
  const StopToken no_stop;
  return generateChild(pMap, no_stop, error) == OperationStatus::success;
}

OperationStatus Node::generateChild(const Polymap* pMap,
                                    const StopToken& stop_token,
                                    std::string* error) {
  C_.clear();
  if (error)
    error->clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (!pMap) {
    if (error)
      *error = "Cannot generate children without a Polymap";
    return OperationStatus::failure;
  }
  if (!seed_is_valid_) {
    if (error)
      *error = "Node seed is outside the representable int range";
    return OperationStatus::failure;
  }
  if (!visibility_region_valid_) {
    if (error)
      *error = visibility_region_error_.empty() ? "Visibility region could not be adapted"
                                                : visibility_region_error_;
    return OperationStatus::failure;
  }
  const Point2d source = {static_cast<double>(seed_.first), static_cast<double>(seed_.second)};
  const auto status = generateChildrenFromSource(pMap,
                                                 Nindex_,
                                                 source,
                                                 exact_geometry::Point(seed_.first, seed_.second),
                                                 start_angle_,
                                                 end_angle_,
                                                 Gcost_,
                                                 visibility_region_,
                                                 C_,
                                                 stop_token,
                                                 error);
  return status;
}

void RaystarCore::outlineMap(std::vector<uint8_t>& costarr, int nx, int ny) {
  const StopToken no_stop;
  (void)outlineMap(costarr, nx, ny, no_stop);
}

OperationStatus RaystarCore::outlineMap(std::vector<uint8_t>& costarr,
                                        int nx,
                                        int ny,
                                        const StopToken& stop_token) {
  if (stop_token.poll())
    return OperationStatus::stopped;
  for (int i = 0; i < nx; i++) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    costarr[i] = 1;
    costarr[(ny - 1) * nx + i] = 1;
  }
  for (int i = 0; i < ny; i++) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    costarr[i * nx] = 1;
    costarr[i * nx + nx - 1] = 1;
  }
  return stop_token.poll() ? OperationStatus::stopped : OperationStatus::success;
}

OperationStatus RaystarCore::getScopedVisibilityRegion(Polymap& theMap,
                                                       Candidate& the_child,
                                                       VisibilityRegion& visibility_region,
                                                       VisibilityRegion& full_visreg,
                                                       const StopToken& stop_token,
                                                       std::string& error) {
  visibility_region.clear();
  full_visreg.clear();
  error.clear();

  const auto fail = [&](const std::string& message) {
    visibility_region.clear();
    full_visreg.clear();
    error = message;
    return OperationStatus::failure;
  };
  const auto stopped = [&]() {
    visibility_region.clear();
    full_visreg.clear();
    error.clear();
    return OperationStatus::stopped;
  };

  if (stop_token.poll())
    return stopped();

  int parent_index = the_child.Nindex_;
  int child_index = the_child.Cindex_;
  if (parent_index < 0 || parent_index >= static_cast<int>(N_.size()))
    return fail("Scoped visibility candidate has an invalid parent index");
  if (child_index < 0 ||
      child_index >= static_cast<int>(N_[static_cast<size_t>(parent_index)].C_.size()))
    return fail("Scoped visibility candidate has an invalid child index");

  const auto& child = N_[static_cast<size_t>(parent_index)].C_[static_cast<size_t>(child_index)];
  const auto* child_vertex = std::get_if<ObstacleVertexId>(&child.c_endpoint_.support);
  if (!child_vertex)
    return fail("Scoped visibility child source is not an exact obstacle vertex");
  std::string endpoint_error;
  auto endpoint_status =
    theMap.validateBoundaryEndpoint(child.c_endpoint_, stop_token, &endpoint_error);
  if (endpoint_status == OperationStatus::stopped)
    return stopped();
  if (endpoint_status == OperationStatus::failure)
    return fail("Scoped visibility child source is invalid: " + endpoint_error);

  const auto new_source_point = child.c_;
  VisibilityRegion fullV;

  std::string visibility_error;
  const auto full_visibility_status = theMap.getVisibilityRegion(
    new_source_point.first, new_source_point.second, fullV, stop_token, &visibility_error);
  if (full_visibility_status == OperationStatus::stopped)
    return stopped();
  if (full_visibility_status == OperationStatus::failure) {
    return fail("Full visibility calculation failed at child source (" +
                std::to_string(new_source_point.first) + ", " +
                std::to_string(new_source_point.second) + "): " + visibility_error);
  }

  BoundaryEndpoint start_obs;
  BoundaryEndpoint end_obs;
  double start_angle, end_angle;
  const std::pair<int, int> child_topology = {child_vertex->obstacle, child_vertex->vertex};

  if (child.is_a_left_gap_) {
    start_obs = child.o_endpoint_;
    auto end_obs_result = theMap.getNextObs(child_topology);
    if (!end_obs_result)
      return fail("Scoped visibility child has an invalid obstacle topology");
    const auto& vertices =
      theMap.obstacles()[static_cast<size_t>(child_vertex->obstacle)].ordered_vertices_;
    const std::pair<int, int> end_topology = {
      child_vertex->obstacle, (child_vertex->vertex + 1) % static_cast<int>(vertices.size())};
    end_obs = exactEndpoint(*end_obs_result, end_topology);
  } else {
    auto start_obs_result = theMap.getPrevObs(child_topology);
    if (!start_obs_result)
      return fail("Scoped visibility child has an invalid obstacle topology");
    const auto& vertices =
      theMap.obstacles()[static_cast<size_t>(child_vertex->obstacle)].ordered_vertices_;
    const std::pair<int, int> start_topology = {
      child_vertex->obstacle,
      (child_vertex->vertex + static_cast<int>(vertices.size()) - 1) %
        static_cast<int>(vertices.size())};
    start_obs = exactEndpoint(*start_obs_result, start_topology);
    end_obs = child.o_endpoint_;
  }
  endpoint_status = theMap.validateBoundaryEndpoint(start_obs, stop_token, &endpoint_error);
  if (endpoint_status == OperationStatus::stopped)
    return stopped();
  if (endpoint_status == OperationStatus::failure)
    return fail("Scoped visibility start boundary is invalid: " + endpoint_error);
  endpoint_status = theMap.validateBoundaryEndpoint(end_obs, stop_token, &endpoint_error);
  if (endpoint_status == OperationStatus::stopped)
    return stopped();
  if (endpoint_status == OperationStatus::failure)
    return fail("Scoped visibility end boundary is invalid: " + endpoint_error);

  start_angle = atan2(start_obs.position.second - new_source_point.second,
                      start_obs.position.first - new_source_point.first);
  end_angle = atan2(end_obs.position.second - new_source_point.second,
                    end_obs.position.first - new_source_point.first);
  const double scope_span = normalize_angle_positive(end_angle - start_angle);
  end_angle = start_angle + scope_span;

  const exact_geometry::Point exact_source(new_source_point.first, new_source_point.second);
  const exact_geometry::Point exact_start = exactPoint(start_obs);
  const exact_geometry::Point exact_end = exactPoint(end_obs);
  if (exact_start == exact_source || exact_end == exact_source)
    return fail("Scoped visibility has a zero-length boundary ray");

  for (size_t i = 0; i < fullV.size(); ++i) {
    if (stop_token.poll())
      return stopped();
    if (exact_geometry::isClosedCounterClockwiseSweepMember(
          exact_source, exact_start, exact_end, exactPoint(fullV[i]))) {
      visibility_region.emplace_back(fullV[i]);
    }
  }

  const auto same_position = [](const BoundaryEndpoint& endpoint,
                                const BoundaryEndpoint& boundary) {
    return exactPoint(endpoint) == exactPoint(boundary);
  };
  size_t loc = 0;
  for (; loc < visibility_region.size(); ++loc) {
    if (stop_token.poll())
      return stopped();
    if (same_position(visibility_region[loc], start_obs))
      break;
  }

  if (loc == visibility_region.size())
    visibility_region.insert(visibility_region.begin(), start_obs);
  else {
    visibility_region.erase(visibility_region.begin(), visibility_region.begin() + loc);
    visibility_region.front() = start_obs;
  }

  loc = 0;
  for (; loc < visibility_region.size(); ++loc) {
    if (stop_token.poll())
      return stopped();
    if (same_position(visibility_region[loc], end_obs))
      break;
  }

  if (loc == visibility_region.size())
    visibility_region.emplace_back(end_obs);
  else {
    visibility_region.erase(visibility_region.begin() + loc + 1, visibility_region.end());
    visibility_region.back() = end_obs;
  }

  std::string validation_error;
  const auto metadata_status =
    theMap.validateVisibilityRegion(visibility_region, stop_token, &validation_error);
  if (metadata_status == OperationStatus::stopped)
    return stopped();
  if (metadata_status == OperationStatus::failure)
    return fail("Scoped " + validation_error);
  const detail::VisibilityGeometryContext geometry_context{
    exact_source, detail::VisibilityBoundaryMode::open_sector, exact_start, exact_end};
  const auto geometry_status = detail::validateVisibilityGeometry(
    visibility_region, geometry_context, stop_token, &validation_error);
  if (geometry_status == OperationStatus::stopped)
    return stopped();
  if (geometry_status == OperationStatus::failure) {
    return fail("Scoped " + validation_error);
  }

  if (stop_token.poll())
    return stopped();
  full_visreg = std::move(fullV);
  return OperationStatus::success;
}

PlanResult RaystarCore::plan(const GridMap& grid_map,
                             int start_x,
                             int start_y,
                             int goal_x,
                             int goal_y,
                             int K,
                             bool allow_self_crossing,
                             const PlanningLimits& limits) {
  std::string legacy_error;
  if (!validatePlanRequest(grid_map, start_x, start_y, goal_x, goal_y, K, limits, legacy_error)) {
    std::vector<Node>().swap(N_);
    PlanResult result;
    result.outcome = PlanningOutcome::invalid_request;
    result.message = std::move(legacy_error);
    return result;
  }
  return plan(
    grid_map,
    PlanEndpoint(start_x, start_y, static_cast<double>(start_x), static_cast<double>(start_y)),
    PlanEndpoint(goal_x, goal_y, static_cast<double>(goal_x), static_cast<double>(goal_y)),
    K,
    allow_self_crossing,
    limits);
}

PlanResult RaystarCore::plan(const GridMap& grid_map,
                             const Point2d& start,
                             const Point2d& goal,
                             int K,
                             bool allow_self_crossing,
                             const PlanningLimits& limits) {
  const auto endpoint_from_point = [&grid_map](const Point2d& point) {
    std::pair<int, int> cell = {-1, -1};
    if (std::isfinite(point.first) && std::isfinite(point.second)) {
      const double x = std::floor(point.first);
      const double y = std::floor(point.second);
      if (x >= static_cast<double>(std::numeric_limits<int>::lowest()) &&
          x <= static_cast<double>(std::numeric_limits<int>::max()) &&
          y >= static_cast<double>(std::numeric_limits<int>::lowest()) &&
          y <= static_cast<double>(std::numeric_limits<int>::max())) {
        cell = {static_cast<int>(x), static_cast<int>(y)};
      }
    }
    (void)grid_map;
    return PlanEndpoint(cell, point);
  };
  return plan(grid_map,
              endpoint_from_point(start),
              endpoint_from_point(goal),
              K,
              allow_self_crossing,
              limits);
}

PlanResult RaystarCore::plan(const GridMap& grid_map,
                             const PlanEndpoint& start,
                             const PlanEndpoint& goal,
                             int K,
                             bool allow_self_crossing,
                             const PlanningLimits& limits) try {
  // Per-plan state must never survive into a new request, including invalid
  // requests and maps for which Polymap reports that no solution exists.
  // Releasing capacity here prevents a previously large request from keeping
  // its entire debug/search allocation resident while the next map is copied.
  std::vector<Node>().swap(N_);

  PlanResult result;
  std::string validation_error;
  if (!validatePlanRequest(grid_map, start, goal, K, limits, validation_error)) {
    result.outcome = PlanningOutcome::invalid_request;
    result.message = std::move(validation_error);
    return result;
  }

  // Every returned path consumes at least the two continuous endpoints.  A
  // bounded reserve avoids repeated geometric growth while still respecting
  // both K and the response-wide path-point budget.
  const size_t path_reserve = std::min(static_cast<size_t>(K), limits.max_path_points / 2);
  result.path_solutions.reserve(path_reserve);

  const int start_x = start.cell_.first;
  const int start_y = start.cell_.second;
  const int goal_x = goal.cell_.first;
  const int goal_y = goal.cell_.second;
  const double start_gx = start.position_.first;
  const double start_gy = start.position_.second;
  const double goal_gx = goal.position_.first;
  const double goal_gy = goal.position_.second;

  using PlanningClock = std::chrono::steady_clock;
  const auto request_start_time = PlanningClock::now();
  const bool timeout_enabled = limits.planning_timeout != std::chrono::milliseconds::max();
  PlanningLimitReached requested_stop_reason = PlanningLimitReached::none;
  const StopToken stop_token([&]() {
    if (limits.cancel_requested && limits.cancel_requested()) {
      requested_stop_reason = PlanningLimitReached::cancelled;
      return true;
    }
    if (timeout_enabled &&
        std::chrono::duration_cast<std::chrono::milliseconds>(
          PlanningClock::now() - request_start_time) >= limits.planning_timeout) {
      requested_stop_reason = PlanningLimitReached::timeout;
      return true;
    }
    return false;
  });

  const auto stop_before_planner = [&](std::shared_ptr<const Polymap> polymap) {
    result.outcome = PlanningOutcome::limit_reached;
    result.limit_reached = requested_stop_reason;
    result.message =
      planningLimitMessage(result.limit_reached, limits, result.path_solutions.size());
    result.polymap = std::move(polymap);
    return result;
  };

  if (stop_token.poll())
    return stop_before_planner(nullptr);

  GridMap work_map = grid_map;
  if (stop_token.poll())
    return stop_before_planner(nullptr);
  if (outlineMap(work_map.data, work_map.width, work_map.height, stop_token) ==
      OperationStatus::stopped) {
    return stop_before_planner(nullptr);
  }

  // Outlining reserves the map border as a geometric boundary.  Reject an
  // endpoint whose cell was consumed by that operation instead of silently
  // clearing it (the former behaviour could turn an occupied start into free
  // space and made a boundary point appear valid).
  if (work_map.at(static_cast<unsigned int>(start_x), static_cast<unsigned int>(start_y)) != 0) {
    result.outcome = PlanningOutcome::invalid_request;
    result.message = "Invalid start: corresponding cell is occupied or on the map boundary";
    return result;
  }
  if (work_map.at(static_cast<unsigned int>(goal_x), static_cast<unsigned int>(goal_y)) != 0) {
    result.outcome = PlanningOutcome::invalid_request;
    result.message = "Invalid goal: corresponding cell is occupied or on the map boundary";
    return result;
  }

  const auto map_start_time = PlanningClock::now();
  auto map_build = Polymap::create(work_map,
                                   start_x,
                                   start_y,
                                   goal_x,
                                   goal_y,
                                   start.position_,
                                   goal.position_,
                                   stop_token,
                                   limits);
  const auto map_end_time = PlanningClock::now();
  result.map_time_ms =
    std::chrono::duration_cast<std::chrono::microseconds>(map_end_time - map_start_time).count() /
    1000.0;

  // Poll after construction as well: a single CGAL primitive cannot be
  // pre-empted, so the deadline may have elapsed inside that call.
  if (map_build.status == PolymapCreateStatus::stopped || stop_token.poll())
    return stop_before_planner(nullptr);

  if (map_build.status == PolymapCreateStatus::no_path) {
    result.success = false;
    result.outcome = PlanningOutcome::no_path;
    result.message = "No path exists between start and goal";
    return result;
  }

  if (map_build.status == PolymapCreateStatus::failure || !map_build) {
    result.success = false;
    // Endpoint validation in the continuous Polymap constructor is an input
    // contract failure, not a generic CDT failure.  Preserve the precise
    // start/goal diagnostic so callers can correct the request directly.
    constexpr const char* invalid_start_prefix = "Invalid start position: ";
    constexpr const char* invalid_goal_prefix = "Invalid goal position: ";
    if (map_build.error.rfind(invalid_start_prefix, 0) == 0) {
      result.outcome = PlanningOutcome::invalid_request;
      result.message = "Invalid start: point must be a strict free-space interior; " +
                       map_build.error.substr(std::string(invalid_start_prefix).size());
    } else if (map_build.error.rfind(invalid_goal_prefix, 0) == 0) {
      result.outcome = PlanningOutcome::invalid_request;
      result.message = "Invalid goal: point must be a strict free-space interior; " +
                       map_build.error.substr(std::string(invalid_goal_prefix).size());
    } else {
      result.outcome = PlanningOutcome::failed;
      result.message = "Map geometry construction failed";
    }
    if (!map_build.error.empty()) {
      if (result.message == "Map geometry construction failed")
        result.message += ": " + map_build.error;
    }
    return result;
  }

  auto theMap = std::make_shared<Polymap>(std::move(*map_build.value));

  const auto planner_start_time = PlanningClock::now();

  const auto fail_planning = [&](const std::string& message,
                                 PlanningOutcome outcome = PlanningOutcome::failed) {
    const auto failure_time = PlanningClock::now();
    result.success = false;
    result.outcome = outcome;
    result.expanded_nodes = N_.size();
    result.message = message;
    result.path_solutions.clear();
    result.plan_time_ms =
      std::chrono::duration_cast<std::chrono::microseconds>(failure_time - planner_start_time)
        .count() /
      1000.0;
    result.polymap = theMap;
    resetSearchState();
    return result;
  };

  std::string endpoint_error;
  const auto start_interior_status =
    theMap->validateFreeSpaceInterior(start.position_, stop_token, &endpoint_error);
  if (start_interior_status == OperationStatus::stopped || stop_token.poll())
    return stop_before_planner(theMap);
  if (start_interior_status == OperationStatus::failure) {
    return fail_planning(
      "Invalid start: point must be a strict free-space interior; " + endpoint_error,
      PlanningOutcome::invalid_request);
  }
  const auto goal_interior_status =
    theMap->validateFreeSpaceInterior(goal.position_, stop_token, &endpoint_error);
  if (goal_interior_status == OperationStatus::stopped || stop_token.poll())
    return stop_before_planner(theMap);
  if (goal_interior_status == OperationStatus::failure) {
    return fail_planning(
      "Invalid goal: point must be a strict free-space interior; " + endpoint_error,
      PlanningOutcome::invalid_request);
  }

  const auto stop_for_limit = [&](PlanningLimitReached limit) {
    const auto stop_time = PlanningClock::now();
    result.outcome = PlanningOutcome::limit_reached;
    result.limit_reached = limit;
    result.expanded_nodes = N_.size();
    result.success = !result.path_solutions.empty();
    result.message = planningLimitMessage(limit, limits, result.path_solutions.size());
    result.plan_time_ms =
      std::chrono::duration_cast<std::chrono::microseconds>(stop_time - planner_start_time)
        .count() /
      1000.0;
    result.polymap = theMap;
    return result;
  };

  const auto stop_for_request = [&]() { return stop_for_limit(requested_stop_reason); };

  size_t accumulated_path_points = 0;

  auto comp = [](const Candidate& a, const Candidate& b) { return a.Fcost_ > b.Fcost_; };

  std::vector<Candidate> queue;
  queue.emplace_back(Candidate(-1, -1, std::hypot(start_gx - goal_gx, start_gy - goal_gy)));

  while (!queue.empty()) {
    if (stop_token.poll())
      return stop_for_request();
    if (N_.size() >= limits.max_nodes)
      return stop_for_limit(PlanningLimitReached::max_nodes);

    Candidate best_candidate = queue[0];
    std::pop_heap(queue.begin(), queue.end(), comp);
    queue.pop_back();
    if (stop_token.poll())
      return stop_for_request();

    int parent_index = best_candidate.Nindex_;
    int child_index = best_candidate.Cindex_;
    std::optional<Node> pending_node;

    if (parent_index == -1) {
      VisibilityRegion Vtemp;
      std::string visibility_error;
      const Point2d start_position = {start_gx, start_gy};
      const auto visibility_status =
        theMap->getRootVisibilityRegion(start_position, Vtemp, stop_token, &visibility_error);
      if (visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (visibility_status == OperationStatus::failure) {
        return fail_planning("Root visibility calculation failed: " + visibility_error);
      }
      if (N_.size() >= limits.max_nodes)
        return stop_for_limit(PlanningLimitReached::max_nodes);

      Node new_node(theMap.get(),
                    0,
                    static_cast<double>(start_x),
                    static_cast<double>(start_y),
                    0.0,
                    best_candidate.Fcost_,
                    Vtemp,
                    stop_token);
      new_node.is_continuous_root_ = true;
      // The root holder remains in N_[0] for the legacy debug/index contract,
      // but its integer seed is not a geometric waypoint.
      new_node.local_shortest_path_.clear();
      if (stop_token.poll())
        return stop_for_request();
      const auto full_visibility_status = new_node.setFullVisibilityRegion(Vtemp, stop_token);
      if (full_visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (full_visibility_status == OperationStatus::failure)
        return fail_planning("Root full visibility projection failed");

      std::string child_error;
      const auto child_status =
        generateChildrenFromSource(theMap.get(),
                                   new_node.Nindex_,
                                   start_position,
                                   exact_geometry::Point(start_gx, start_gy),
                                   new_node.start_angle_,
                                   new_node.end_angle_,
                                   new_node.Gcost_,
                                   new_node.visibility_region_,
                                   new_node.C_,
                                   stop_token,
                                   &child_error);
      if (child_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (child_status == OperationStatus::failure)
        return fail_planning("Root child generation failed: " + child_error);
      if (N_.size() >= limits.max_nodes)
        return stop_for_limit(PlanningLimitReached::max_nodes);
      pending_node.emplace(std::move(new_node));
    } else {
      if (parent_index < 0 || parent_index >= static_cast<int>(N_.size()))
        return fail_planning("Planner candidate has an invalid parent index");
      if (child_index < 0 ||
          child_index >= static_cast<int>(N_[static_cast<size_t>(parent_index)].C_.size()))
        return fail_planning("Planner candidate has an invalid child index");

      auto new_source_point = N_[parent_index].C_[child_index].c_;
      int new_node_index = static_cast<int>(N_.size());

      bool path_self_crosses = false;
      if (!allow_self_crossing) {
        const auto crossing_status =
          hybridPathSelfCrosses({start_gx, start_gy},
                                N_[parent_index].local_shortest_path_,
                                {static_cast<double>(new_source_point.first),
                                 static_cast<double>(new_source_point.second)},
                                path_self_crosses,
                                stop_token);
        if (crossing_status == OperationStatus::stopped || stop_token.poll())
          return stop_for_request();
        if (crossing_status == OperationStatus::failure)
          return fail_planning("Path self-crossing validation failed");
      }
      if (path_self_crosses)
        continue;

      VisibilityRegion Vtemp;
      VisibilityRegion fullVtemp;
      std::string visibility_error;
      const auto visibility_status = getScopedVisibilityRegion(
        *theMap, best_candidate, Vtemp, fullVtemp, stop_token, visibility_error);
      if (visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (visibility_status == OperationStatus::failure) {
        return fail_planning("Scoped visibility calculation failed: " + visibility_error);
      }
      if (N_.size() >= limits.max_nodes)
        return stop_for_limit(PlanningLimitReached::max_nodes);

      Node new_node(theMap.get(),
                    new_node_index,
                    new_source_point.first,
                    new_source_point.second,
                    N_[parent_index].C_[child_index].c_gcost_,
                    N_[parent_index].C_[child_index].c_hcost_,
                    parent_index,
                    Vtemp,
                    stop_token);
      if (stop_token.poll())
        return stop_for_request();
      const auto full_visibility_status = new_node.setFullVisibilityRegion(fullVtemp, stop_token);
      if (full_visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (full_visibility_status == OperationStatus::failure)
        return fail_planning("Scoped full visibility projection failed");

      new_node.local_shortest_path_.clear();
      new_node.local_shortest_path_.reserve(N_[parent_index].local_shortest_path_.size() + 1);
      for (const auto& waypoint : N_[parent_index].local_shortest_path_) {
        if (stop_token.poll())
          return stop_for_request();
        new_node.local_shortest_path_.emplace_back(waypoint);
      }
      new_node.local_shortest_path_.emplace_back(new_source_point);
      new_node.path_node_index_ = N_[parent_index].path_node_index_;
      new_node.path_node_index_.emplace_back(new_node_index);
      new_node.start_angle_ = N_[parent_index].C_[child_index].start_angle_;
      new_node.end_angle_ = N_[parent_index].C_[child_index].end_angle_;

      std::string child_error;
      const auto child_status = new_node.generateChild(theMap.get(), stop_token, &child_error);
      if (child_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (child_status == OperationStatus::failure)
        return fail_planning("Scoped child generation failed: " + child_error);
      if (N_.size() >= limits.max_nodes)
        return stop_for_limit(PlanningLimitReached::max_nodes);
      pending_node.emplace(std::move(new_node));
    }

    if (!pending_node)
      return fail_planning("Planner failed to construct the selected node");

    for (auto& child : pending_node->C_) {
      if (stop_token.poll())
        return stop_for_request();
      child.c_hcost_ = std::hypot(child.c_endpoint_.position.first - goal_gx,
                                  child.c_endpoint_.position.second - goal_gy);
    }
    if (stop_token.poll())
      return stop_for_request();

    N_.emplace_back(std::move(*pending_node));
    for (const auto& child : N_.back().C_) {
      if (stop_token.poll())
        return stop_for_request();
      queue.emplace_back(Candidate(child.Nindex_, child.Cindex_, child.c_gcost_ + child.c_hcost_));
      std::push_heap(queue.begin(), queue.end(), comp);
    }

    if (stop_token.poll())
      return stop_for_request();
    const auto point_location_mode = visibilityBoundaryModeForNode(*theMap, N_.back());
    const exact_geometry::Point exact_node_source =
      N_.back().is_continuous_root_
        ? exact_geometry::Point(start_gx, start_gy)
        : exact_geometry::Point(N_.back().seed_.first, N_.back().seed_.second);
    std::pair<bool, bool> b = {false, false};
    const auto point_location_status =
      detail::classifyPointInVisibilityBoundary(N_.back().visibility_region_,
                                                exact_geometry::Point(goal_gx, goal_gy),
                                                exact_node_source,
                                                point_location_mode,
                                                b,
                                                stop_token);
    if (point_location_status == OperationStatus::stopped || stop_token.poll())
      return stop_for_request();
    if (point_location_status == OperationStatus::failure)
      return fail_planning("Goal point-location validation failed");
    if (b.first || b.second) {
      bool goal_segment_crosses = false;
      if (!allow_self_crossing) {
        const auto crossing_status = hybridPathSelfCrosses({start_gx, start_gy},
                                                           N_.back().local_shortest_path_,
                                                           {goal_gx, goal_gy},
                                                           goal_segment_crosses,
                                                           stop_token);
        if (crossing_status == OperationStatus::stopped || stop_token.poll())
          return stop_for_request();
        if (crossing_status == OperationStatus::failure)
          return fail_planning("Goal path self-crossing validation failed");
      }

      if (!goal_segment_crosses) {
        const size_t turning_point_count = N_.back().local_shortest_path_.size();
        // PathSolution always contains the two continuous endpoints in
        // addition to its integer turning points.  Check both the per-path
        // addition and the response-wide total before constructing/copying a
        // solution, so an adversarial search cannot grow an unbounded result.
        if (turning_point_count > limits.max_path_points - 2 ||
            accumulated_path_points > limits.max_path_points - (turning_point_count + 2)) {
          return stop_for_limit(PlanningLimitReached::max_path_points);
        }
        const Point2d node_source = N_.back().is_continuous_root_
                                      ? Point2d{start_gx, start_gy}
                                      : Point2d{static_cast<double>(N_.back().seed_.first),
                                                static_cast<double>(N_.back().seed_.second)};
        const double new_path_length =
          N_.back().Gcost_ + std::hypot(node_source.first - goal_gx, node_source.second - goal_gy);

        PathSolution solution({start_gx, start_gy},
                              N_.back().local_shortest_path_,
                              {goal_gx, goal_gy},
                              new_path_length,
                              N_.back().path_node_index_);
        if (stop_token.poll())
          return stop_for_request();
        result.path_solutions.emplace_back(std::move(solution));
        accumulated_path_points += turning_point_count + 2;
        if (stop_token.poll())
          return stop_for_request();
        if (static_cast<int>(result.path_solutions.size()) >= K)
          break;
      }
    }
  }

  if (stop_token.poll())
    return stop_for_request();

  const auto planner_end_time = PlanningClock::now();
  result.plan_time_ms =
    std::chrono::duration_cast<std::chrono::microseconds>(planner_end_time - planner_start_time)
      .count() /
    1000.0;

  result.success = !result.path_solutions.empty();
  result.expanded_nodes = N_.size();
  if (!result.success) {
    result.outcome = PlanningOutcome::no_path;
    result.message = "Planning completed but no path found";
  } else {
    result.outcome = PlanningOutcome::complete;
  }
  result.polymap = theMap;
  return result;
} catch (const std::exception& exception) {
  resetSearchState();
  PlanResult result;
  result.message = boundedPlanningExceptionMessage(exception.what());
  return result;
} catch (...) {
  resetSearchState();
  PlanResult result;
  result.message = "Planning failed with an unknown exception";
  return result;
}

}  // namespace raystar
