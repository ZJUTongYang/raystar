#include <raystar/raystar_core.h>
#include <chrono>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <limits>
#include <exception>
#include <optional>

#include "exact_point_location.h"
#include "conservative_path_length.h"
#include "visibility_validation.h"

#include "raystar_core_detail.h"

namespace raystar {

using namespace core_impl;

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


}  // namespace raystar
