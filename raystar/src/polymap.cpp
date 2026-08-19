#include <raystar/polymap.h>
#include "visibility_validation.h"
#include <unordered_set>
#include <stack>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <functional>
#include <list>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <CGAL/exceptions.h>
#include <CGAL/number_utils.h>

#include "polymap_detail.h"

namespace raystar {

using namespace polymap_impl;

bool Polymap::getPolyObstacles(int start_x, int start_y, int goal_x, int goal_y) {
  const StopToken no_stop;
  return getPolyObstacles(start_x, start_y, goal_x, goal_y, no_stop) == OperationStatus::success;
}

OperationStatus Polymap::getPolyObstacles(
  int start_x, int start_y, int goal_x, int goal_y, const StopToken& stop_token) {
  return getPolyObstacles(
    start_x,
    start_y,
    std::vector<PolymapEndpoint>{
      {goal_x, goal_y, {static_cast<double>(goal_x), static_cast<double>(goal_y)}}},
    stop_token);
}

OperationStatus Polymap::getPolyObstacles(int start_x,
                                          int start_y,
                                          const std::vector<PolymapEndpoint>& goals,
                                          const StopToken& stop_token) {
  return getPolyObstacles(start_x, start_y, goals, stop_token, std::nullopt);
}

OperationStatus Polymap::getPolyObstacles(
  int start_x,
  int start_y,
  const std::vector<PolymapEndpoint>& goals,
  const StopToken& stop_token,
  std::optional<size_t> max_raw_contour_vertices) {
  if (stop_token.poll())
    return OperationStatus::stopped;
  construction_error_.clear();
  if (getPolyObstaclesImpl(
        start_x, start_y, goals, stop_token, max_raw_contour_vertices)) {
    // obs_ now contains a newly extracted, unsimplified contour set.  Any
    // topology registry, CDT, facet table or visibility cache from a previous
    // build refers to different vertex indices and must not remain usable.
    solution_exist_ = true;
    construction_stopped_ = false;
    construction_error_.clear();
    std::fill(vertices_location_x_flat_.begin(), vertices_location_x_flat_.end(), -1);
    std::fill(vertices_location_y_flat_.begin(), vertices_location_y_flat_.end(), -1);
    clearCGALRelatedState();
    return OperationStatus::success;
  }
  if (stop_token.poll())
    return OperationStatus::stopped;
  // A normal refresh failure invalidates every object derived from the old
  // obstacle set; keeping a ready CDT/cache next to an empty obs_ vector would
  // expose stale geometry.  A cooperative stop is different: the previously
  // committed state remains untouched.
  solution_exist_ = false;
  construction_stopped_ = false;
  obs_.clear();
  std::fill(vertices_location_x_flat_.begin(), vertices_location_x_flat_.end(), -1);
  std::fill(vertices_location_y_flat_.begin(), vertices_location_y_flat_.end(), -1);
  clearCGALRelatedState();
  return OperationStatus::failure;
}

bool Polymap::getPolyObstaclesImpl(
  int start_x,
  int start_y,
  int goal_x,
  int goal_y,
  const StopToken& stop_token,
  std::optional<size_t> max_raw_contour_vertices) {
  return getPolyObstaclesImpl(
    start_x,
    start_y,
    std::vector<PolymapEndpoint>{
      {goal_x, goal_y, {static_cast<double>(goal_x), static_cast<double>(goal_y)}}},
    stop_token,
    max_raw_contour_vertices);
}

bool Polymap::getPolyObstaclesImpl(int start_x,
                                   int start_y,
                                   const std::vector<PolymapEndpoint>& goals,
                                   const StopToken& stop_token,
                                   std::optional<size_t> max_raw_contour_vertices) {
  if (stop_token.poll())
    return false;
  if (xsize_ <= 0 || ysize_ <= 0)
    return false;
  const size_t width = static_cast<size_t>(xsize_);
  const size_t height = static_cast<size_t>(ysize_);
  if (width > std::numeric_limits<size_t>::max() / height)
    return false;
  const size_t cell_count = width * height;
  if (cell_count > static_cast<size_t>(std::numeric_limits<int>::max()) ||
      data_.size() != cell_count) {
    return false;
  }
  const auto in_bounds = [this](int x, int y) {
    return x >= 0 && y >= 0 && x < xsize_ && y < ysize_;
  };
  if (!in_bounds(start_x, start_y) || goals.empty())
    return false;
  for (const auto& goal : goals) {
    if (!in_bounds(goal.cell_x, goal.cell_y))
      return false;
  }

  std::vector<Obs> candidate_obstacles;
  unsigned int nx = xsize_;
  unsigned int ny = ysize_;

  std::vector<int> mask(nx * ny, 0);

  std::unordered_map<UndirectedGridEdgeKey, DirectedGridEdge, UndirectedGridEdgeKeyHash> edges;
  const auto add_directed_edge = [&](int from, int to) {
    const auto insertion =
      edges.emplace(makeUndirectedGridEdgeKey(from, to), DirectedGridEdge{from, to});
    if (max_raw_contour_vertices && edges.size() > *max_raw_contour_vertices) {
      if (insertion.second)
        edges.erase(insertion.first);
      construction_error_ =
        "Reference-shortening raw contour exceeds the vertex budget " +
        std::to_string(*max_raw_contour_vertices);
      return false;
    }
    return true;
  };
  std::stack<int> Q;
  Q.push(start_x + start_y * nx);
  while (!Q.empty()) {
    if (stop_token.poll())
      return false;
    int cur = Q.top();
    Q.pop();
    int x = cur % nx;
    int y = (cur - x) / nx;
    if (data_[cur] != 0 || mask[x + y * nx] != 0)
      continue;

    if (cur - 1 >= 0 && data_[cur - 1] != 0 && !add_directed_edge(cur, cur + nx))
      return false;
    if (cur + 1 < static_cast<int>(nx * ny) && data_[cur + 1] != 0 &&
        !add_directed_edge(cur + nx + 1, cur + 1)) {
      return false;
    }
    if (cur - static_cast<int>(nx) >= 0 && data_[cur - nx] != 0 &&
        !add_directed_edge(cur + 1, cur)) {
      return false;
    }
    if (cur + static_cast<int>(nx) < static_cast<int>(nx * ny) && data_[cur + nx] != 0 &&
        !add_directed_edge(cur + nx, cur + nx + 1)) {
      return false;
    }

    mask[x + y * nx] = 1;
    if (x > 0 && data_[cur - 1] == 0)
      Q.push(cur - 1);
    if (x < static_cast<int>(nx) - 1 && data_[cur + 1] == 0)
      Q.push(cur + 1);
    if (y > 0 && data_[cur - nx] == 0)
      Q.push(cur - nx);
    if (y < static_cast<int>(ny) - 1 && data_[cur + nx] == 0)
      Q.push(cur + nx);
  }

  for (const auto& goal : goals) {
    if (mask[static_cast<size_t>(goal.cell_x) + static_cast<size_t>(goal.cell_y) * nx] == 0)
      return false;
  }

  std::list<std::pair<int, int>> boundary_points;
  while (!edges.empty()) {
    if (stop_token.poll())
      return false;
    candidate_obstacles.emplace_back(Obs());
    boundary_points.clear();

    auto first_iter = edges.begin();
    const int prev = first_iter->second.from;
    const int next = first_iter->second.to;
    int prev_x = prev % nx;
    int prev_y = (prev - prev_x) / nx;
    int next_x = next % nx;
    int next_y = (next - next_x) / nx;

    boundary_points.emplace_back(prev_x, prev_y);
    boundary_points.emplace_back(next_x, next_y);

    auto biter = boundary_points.end();
    int x, y;
    while (1) {
      if (stop_token.poll())
        return false;
      x = boundary_points.back().first;
      y = boundary_points.back().second;
      int cur = x + y * nx;

      int lb_free = 0, lt_free = 0, rb_free = 0, rt_free = 0;
      if (x > 0 && y > 0 && data_[cur - nx - 1] == 0)
        lb_free = 1;
      if (x > 0 && data_[cur - 1] == 0)
        lt_free = 1;
      if (y > 0 && data_[cur - nx] == 0)
        rb_free = 1;
      if (data_[cur] == 0)
        rt_free = 1;

      int num = lb_free * 8 + lt_free * 4 + rb_free * 2 + rt_free;

      switch (num) {
      case 1:
        boundary_points.emplace_back(x, y + 1);
        break;
      case 2:
        boundary_points.emplace_back(x + 1, y);
        break;
      case 3:
        boundary_points.emplace_back(x, y + 1);
        break;
      case 4:
        boundary_points.emplace_back(x - 1, y);
        break;
      case 5:
        boundary_points.emplace_back(x - 1, y);
        break;
      case 6:
        biter = std::prev(boundary_points.end(), 2);
        if ((*biter).first == x && (*biter).second == y - 1)
          boundary_points.emplace_back(x + 1, y);
        else if ((*biter).first == x && (*biter).second == y + 1)
          boundary_points.emplace_back(x - 1, y);
        break;
      case 7:
        boundary_points.emplace_back(x - 1, y);
        break;
      case 8:
        boundary_points.emplace_back(x, y - 1);
        break;
      case 9:
        biter = std::prev(boundary_points.end(), 2);
        if ((*biter).first == x + 1 && (*biter).second == y)
          boundary_points.emplace_back(x, y + 1);
        else if ((*biter).first == x - 1 && (*biter).second == y)
          boundary_points.emplace_back(x, y - 1);
        break;
      case 10:
        boundary_points.emplace_back(x + 1, y);
        break;
      case 11:
        boundary_points.emplace_back(x, y + 1);
        break;
      case 12:
        boundary_points.emplace_back(x, y - 1);
        break;
      case 13:
        boundary_points.emplace_back(x, y - 1);
        break;
      case 14:
        boundary_points.emplace_back(x + 1, y);
        break;
      }

      if (boundary_points.back().first == boundary_points.front().first &&
          boundary_points.back().second == boundary_points.front().second) {
        boundary_points.pop_back();
        break;
      }
    }

    for (auto& bp : boundary_points) {
      if (stop_token.poll())
        return false;
      candidate_obstacles.back().ordered_vertices_.emplace_back(bp);
    }

    for (auto it = candidate_obstacles.back().ordered_vertices_.begin();
         it != candidate_obstacles.back().ordered_vertices_.end();
         ++it) {
      if (stop_token.poll())
        return false;
      int curr = it->first + it->second * nx;
      int nxt;
      if (it == candidate_obstacles.back().ordered_vertices_.end() - 1)
        nxt = candidate_obstacles.back().ordered_vertices_.front().first +
              candidate_obstacles.back().ordered_vertices_.front().second * nx;
      else
        nxt = std::next(it)->first + std::next(it)->second * nx;
      edges.erase(makeUndirectedGridEdgeKey(curr, nxt));
    }
  }
  if (stop_token.poll())
    return false;
  obs_ = std::move(candidate_obstacles);
  return true;
}

std::pair<int, int> Polymap::locateVertex(int x, int y) const {
  if (xsize_ <= 0 || ysize_ <= 0 || x < 0 || y < 0 || x >= xsize_ || y >= ysize_)
    return {-1, -1};

  const size_t idx = static_cast<size_t>(y) * static_cast<size_t>(xsize_) + static_cast<size_t>(x);
  if (idx >= vertices_location_x_flat_.size() || idx >= vertices_location_y_flat_.size())
    return {-1, -1};

  return {vertices_location_x_flat_[idx], vertices_location_y_flat_[idx]};
}

std::pair<int, int> Polymap::locateVertex(double x, double y) const {
  if (!std::isfinite(x) || !std::isfinite(y))
    return {-1, -1};

  const double rounded_x = std::round(x);
  const double rounded_y = std::round(y);
  constexpr double integer_tolerance = 1e-9;
  if (std::abs(x - rounded_x) > integer_tolerance || std::abs(y - rounded_y) > integer_tolerance ||
      rounded_x < 0.0 || rounded_y < 0.0 || rounded_x >= static_cast<double>(xsize_) ||
      rounded_y >= static_cast<double>(ysize_)) {
    return {-1, -1};
  }

  return locateVertex(static_cast<int>(rounded_x), static_cast<int>(rounded_y));
}

bool Polymap::validateObstacleTopology(std::string& error, bool validate_edge_relations) const {
  const StopToken no_stop;
  return validateObstacleTopologyImpl(error, validate_edge_relations, no_stop);
}

OperationStatus Polymap::validateObstacleTopology(std::string& error,
                                                  bool validate_edge_relations,
                                                  const StopToken& stop_token) const {
  error.clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (validateObstacleTopologyImpl(error, validate_edge_relations, stop_token)) {
    if (stop_token.poll()) {
      error.clear();
      return OperationStatus::stopped;
    }
    return OperationStatus::success;
  }
  if (stop_token.poll()) {
    error.clear();
    return OperationStatus::stopped;
  }
  return OperationStatus::failure;
}

bool Polymap::validateObstacleTopologyImpl(std::string& error,
                                           bool validate_edge_relations,
                                           const StopToken& stop_token) const {
  error.clear();
  const auto fail = [&](std::string message) {
    error = "Unsupported obstacle topology: " + std::move(message);
    return false;
  };

  if (stop_token.poll())
    return false;

  struct ObstacleEdgeRef {
    size_t obstacle;
    size_t edge;
    size_t obstacle_size;
    std::pair<int, int> from_coordinate;
    std::pair<int, int> to_coordinate;
    exact_geometry::Point from;
    exact_geometry::Point to;
  };

  std::unordered_map<uint64_t, ObstacleVertexId> vertex_owners;
  std::vector<ObstacleEdgeRef> edges;
  for (size_t obstacle_index = 0; obstacle_index < obs_.size(); ++obstacle_index) {
    if (stop_token.poll())
      return false;
    const auto& vertices = obs_[obstacle_index].ordered_vertices_;
    if (vertices.size() < 3) {
      return fail("obstacle " + std::to_string(obstacle_index) + " has fewer than three vertices");
    }

    exact_geometry::FT twice_signed_area(0);
    for (size_t vertex_index = 0; vertex_index < vertices.size(); ++vertex_index) {
      if (stop_token.poll())
        return false;
      const auto& vertex = vertices[vertex_index];
      const auto& next = vertices[(vertex_index + 1) % vertices.size()];
      if (vertex.first < 0 || vertex.second < 0 || vertex.first >= xsize_ ||
          vertex.second >= ysize_) {
        return fail("obstacle " + std::to_string(obstacle_index) + " vertex " +
                    std::to_string(vertex_index) + " at (" + std::to_string(vertex.first) + ", " +
                    std::to_string(vertex.second) + ") is outside the vertex registry");
      }

      const uint64_t vertex_key =
        (static_cast<uint64_t>(static_cast<uint32_t>(vertex.first)) << 32) |
        static_cast<uint32_t>(vertex.second);
      const ObstacleVertexId current_owner{static_cast<int>(obstacle_index),
                                           static_cast<int>(vertex_index)};
      const auto inserted = vertex_owners.emplace(vertex_key, current_owner);
      if (!inserted.second) {
        const auto& previous_owner = inserted.first->second;
        if (previous_owner.obstacle != current_owner.obstacle) {
          return fail("obstacles " + std::to_string(previous_owner.obstacle) + " and " +
                      std::to_string(current_owner.obstacle) + " share vertex (" +
                      std::to_string(vertex.first) + ", " + std::to_string(vertex.second) +
                      "); zero-clearance obstacle contacts are not supported");
        }
        return fail("obstacle " + std::to_string(current_owner.obstacle) + " repeats vertex (" +
                    std::to_string(vertex.first) + ", " + std::to_string(vertex.second) +
                    ") at indices " + std::to_string(previous_owner.vertex) + " and " +
                    std::to_string(current_owner.vertex) +
                    "; self-touching contours are not supported");
      }

      if (vertex == next) {
        return fail("obstacle " + std::to_string(obstacle_index) + " edge " +
                    std::to_string(vertex_index) + " has zero length at (" +
                    std::to_string(vertex.first) + ", " + std::to_string(vertex.second) + ")");
      }

      twice_signed_area += exact_geometry::FT(vertex.first) * exact_geometry::FT(next.second) -
                           exact_geometry::FT(next.first) * exact_geometry::FT(vertex.second);
      if (validate_edge_relations) {
        edges.push_back(ObstacleEdgeRef{obstacle_index,
                                        vertex_index,
                                        vertices.size(),
                                        vertex,
                                        next,
                                        exact_geometry::Point(vertex.first, vertex.second),
                                        exact_geometry::Point(next.first, next.second)});
      }
    }

    if (twice_signed_area == exact_geometry::FT(0)) {
      return fail("obstacle " + std::to_string(obstacle_index) + " has zero signed area");
    }
  }

  if (!validate_edge_relations)
    return true;

  for (size_t first_index = 0; first_index < edges.size(); ++first_index) {
    if (stop_token.poll())
      return false;
    const auto& first = edges[first_index];
    for (size_t second_index = first_index + 1; second_index < edges.size(); ++second_index) {
      if (stop_token.poll())
        return false;
      const auto& second = edges[second_index];
      const auto relation = classifyExactSegments(first.from, first.to, second.from, second.to);

      const bool same_obstacle = first.obstacle == second.obstacle;
      const bool first_precedes_second =
        same_obstacle && (first.edge + 1) % first.obstacle_size == second.edge;
      const bool second_precedes_first =
        same_obstacle && (second.edge + 1) % second.obstacle_size == first.edge;
      if (first_precedes_second || second_precedes_first) {
        const exact_geometry::Point& expected_contact =
          first_precedes_second ? first.to : first.from;
        if (relation.relation == ExactSegmentRelation::endpoint_touch && relation.contact &&
            *relation.contact == expected_contact) {
          continue;
        }
      } else if (relation.relation == ExactSegmentRelation::disjoint) {
        continue;
      }

      std::ostringstream message;
      message << "obstacle " << first.obstacle << " edge " << first.edge << " [("
              << first.from_coordinate.first << ", " << first.from_coordinate.second << ") -> ("
              << first.to_coordinate.first << ", " << first.to_coordinate.second
              << ")] and obstacle " << second.obstacle << " edge " << second.edge << " [("
              << second.from_coordinate.first << ", " << second.from_coordinate.second << ") -> ("
              << second.to_coordinate.first << ", " << second.to_coordinate.second << ")] "
              << segmentRelationDescription(relation.relation);
      return fail(message.str());
    }
  }

  return true;
}

bool Polymap::simplificationChordIsTopologicallySafe(size_t obstacle_index,
                                                     size_t current_index) const {
  const StopToken no_stop;
  return simplificationChordIsTopologicallySafeImpl(obstacle_index, current_index, no_stop);
}

OperationStatus Polymap::simplificationChordIsTopologicallySafe(size_t obstacle_index,
                                                                size_t current_index,
                                                                bool& safe,
                                                                const StopToken& stop_token) const {
  safe = false;
  if (stop_token.poll())
    return OperationStatus::stopped;
  safe = simplificationChordIsTopologicallySafeImpl(obstacle_index, current_index, stop_token);
  if (stop_token.poll()) {
    safe = false;
    return OperationStatus::stopped;
  }
  return OperationStatus::success;
}

bool Polymap::simplificationChordIsTopologicallySafeImpl(size_t obstacle_index,
                                                         size_t current_index,
                                                         const StopToken& stop_token) const {
  if (stop_token.poll())
    return false;
  if (obstacle_index >= obs_.size())
    return false;
  const auto& vertices = obs_[obstacle_index].ordered_vertices_;
  if (vertices.size() <= 3 || current_index >= vertices.size())
    return false;

  const size_t previous_index = (current_index + vertices.size() - 1) % vertices.size();
  const size_t next_index = (current_index + 1) % vertices.size();
  const auto& previous_coordinate = vertices[previous_index];
  const auto& next_coordinate = vertices[next_index];
  if (previous_coordinate == next_coordinate)
    return false;

  const exact_geometry::Point previous(previous_coordinate.first, previous_coordinate.second);
  const exact_geometry::Point next(next_coordinate.first, next_coordinate.second);
  const size_t edge_before_previous = (previous_index + vertices.size() - 1) % vertices.size();

  for (size_t other_obstacle = 0; other_obstacle < obs_.size(); ++other_obstacle) {
    if (stop_token.poll())
      return false;
    const auto& other_vertices = obs_[other_obstacle].ordered_vertices_;
    for (size_t edge_index = 0; edge_index < other_vertices.size(); ++edge_index) {
      if (stop_token.poll())
        return false;
      if (other_obstacle == obstacle_index &&
          (edge_index == previous_index || edge_index == current_index)) {
        continue;
      }

      const auto& edge_from_coordinate = other_vertices[edge_index];
      const auto& edge_to_coordinate = other_vertices[(edge_index + 1) % other_vertices.size()];
      if (edge_from_coordinate == edge_to_coordinate)
        return false;
      const exact_geometry::Point edge_from(edge_from_coordinate.first,
                                            edge_from_coordinate.second);
      const exact_geometry::Point edge_to(edge_to_coordinate.first, edge_to_coordinate.second);
      const auto relation = classifyExactSegments(previous, next, edge_from, edge_to);
      if (relation.relation == ExactSegmentRelation::disjoint)
        continue;

      std::optional<exact_geometry::Point> allowed_contact;
      if (other_obstacle == obstacle_index) {
        if (edge_index == edge_before_previous)
          allowed_contact = previous;
        else if (edge_index == next_index)
          allowed_contact = next;
      }
      if (allowed_contact && relation.relation == ExactSegmentRelation::endpoint_touch &&
          relation.contact && *relation.contact == *allowed_contact) {
        continue;
      }
      return false;
    }
  }
  return true;
}

PolymapCreateResult Polymap::create(const GridMap& grid_map,
                                    int start_x,
                                    int start_y,
                                    int goal_x,
                                    int goal_y,
                                    const PlanningLimits& limits) {
  return create(grid_map, start_x, start_y, goal_x, goal_y, StopToken{}, limits);
}

PolymapCreateResult Polymap::create(const GridMap& grid_map,
                                    int start_x,
                                    int start_y,
                                    int goal_x,
                                    int goal_y,
                                    const StopToken& stop_token,
                                    const PlanningLimits& limits) {
  PolymapCreateResult result;
  MapResourceEstimate estimate;
  if (!validateMapResourceBudget(static_cast<size_t>(grid_map.width),
                                 static_cast<size_t>(grid_map.height),
                                 grid_map.data.size(),
                                 limits,
                                 estimate,
                                 result.error)) {
    return result;
  }
  (void)estimate;
  return finishCreation(Polymap(grid_map, start_x, start_y, goal_x, goal_y, stop_token));
}

PolymapCreateResult Polymap::finishCreation(Polymap&& candidate) {
  PolymapCreateResult result;
  if (candidate.construction_stopped_) {
    result.status = PolymapCreateStatus::stopped;
    return result;
  }
  if (candidate.no_path_) {
    result.status = PolymapCreateStatus::no_path;
    result.error = std::move(candidate.construction_error_);
    return result;
  }
  if (!candidate.solution_exist_ || !candidate.cdt_ready_) {
    result.status = PolymapCreateStatus::failure;
    result.error = candidate.construction_error_.empty() ? "Map geometry construction failed"
                                                         : std::move(candidate.construction_error_);
    return result;
  }

  result.status = PolymapCreateStatus::ready;
  result.value.emplace(std::move(candidate));
  return result;
}

Polymap::Polymap(const GridMap& grid_map,
                 int start_x,
                 int start_y,
                 int goal_x,
                 int goal_y,
                 const StopToken& stop_token)
  : xsize_(0), ysize_(0) {
  const auto stop_construction = [&]() { clearStoppedConstructionState(); };

  if (stop_token.poll()) {
    stop_construction();
    return;
  }

  size_t cell_count = 0;
  if (!prepareGridStorage(grid_map, xsize_, ysize_, cell_count, construction_error_)) {
    return;
  }
  data_ = grid_map.data;
  vertices_location_x_flat_.resize(cell_count, -1);
  vertices_location_y_flat_.resize(cell_count, -1);

  if (start_x < 0 || start_y < 0 || goal_x < 0 || goal_y < 0 || start_x >= xsize_ ||
      start_y >= ysize_ || goal_x >= xsize_ || goal_y >= ysize_) {
    construction_error_ = "Start or goal grid cell is outside map bounds";
    return;
  }
  const size_t start_index =
    static_cast<size_t>(start_x) + static_cast<size_t>(start_y) * static_cast<size_t>(xsize_);
  const size_t goal_index =
    static_cast<size_t>(goal_x) + static_cast<size_t>(goal_y) * static_cast<size_t>(xsize_);
  if (data_[start_index] != 0 || data_[goal_index] != 0) {
    construction_error_ = "Start or goal grid cell is occupied";
    return;
  }

  const OperationStatus obstacle_status =
    getPolyObstacles(start_x, start_y, goal_x, goal_y, stop_token);
  if (obstacle_status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  solution_exist_ = obstacle_status == OperationStatus::success;
  if (!solution_exist_) {
    no_path_ = true;
    return;
  }

  OperationStatus status = OperationStatus::success;
  const bool has_outer_contour = hasClockwiseOuterContour(obs_, stop_token);
  if (stop_token.poll()) {
    stop_construction();
    return;
  }

  // Validate the caller's integer endpoints against the raw reachable
  // contours before any vertex simplification can move a boundary edge.  The
  // same check is repeated after simplification below because the simplified
  // contour is the geometry consumed by CDT/visibility.
  std::string endpoint_error;
  if (has_outer_contour) {
    status =
      validateFreeSpaceInterior(Point2d{static_cast<double>(start_x), static_cast<double>(start_y)},
                                stop_token,
                                &endpoint_error);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid start position: " + endpoint_error;
      return;
    }
    status =
      validateFreeSpaceInterior(Point2d{static_cast<double>(goal_x), static_cast<double>(goal_y)},
                                stop_token,
                                &endpoint_error);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid goal position: " + endpoint_error;
      return;
    }
  }

  // Raw contours consist only of unique unit grid edges. An O(V) incidence
  // pass is sufficient here to reject checkerboard/saddle contacts before
  // simplification, without paying an unnecessary O(E^2) exact edge scan.
  status = validateObstacleTopology(construction_error_, false, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;

  status = simplifyPolyObstacles(start_x, start_y, goal_x, goal_y, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure) {
    construction_error_ = "Obstacle simplification failed";
    return;
  }

  status = validateObstacleTopology(construction_error_, true, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;

  if (has_outer_contour) {
    status =
      validateFreeSpaceInterior(Point2d{static_cast<double>(start_x), static_cast<double>(start_y)},
                                stop_token,
                                &endpoint_error);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid start position: " + endpoint_error;
      return;
    }
    status =
      validateFreeSpaceInterior(Point2d{static_cast<double>(goal_x), static_cast<double>(goal_y)},
                                stop_token,
                                &endpoint_error);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid goal position: " + endpoint_error;
      return;
    }
  }

  status = registerVertices(construction_error_, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;

  status = constructCGALRelated(&validateCDTWithCGAL, construction_error_, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;
}

PolymapCreateResult Polymap::create(const GridMap& grid_map,
                                    int start_x,
                                    int start_y,
                                    int goal_x,
                                    int goal_y,
                                    const Point2d& start_position,
                                    const Point2d& goal_position,
                                    const PlanningLimits& limits) {
  return create(
    grid_map, start_x, start_y, goal_x, goal_y, start_position, goal_position, StopToken{}, limits);
}

PolymapCreateResult Polymap::create(const GridMap& grid_map,
                                    int start_x,
                                    int start_y,
                                    int goal_x,
                                    int goal_y,
                                    const Point2d& start_position,
                                    const Point2d& goal_position,
                                    const StopToken& stop_token,
                                    const PlanningLimits& limits) {
  PolymapCreateResult result;
  MapResourceEstimate estimate;
  if (!validateMapResourceBudget(static_cast<size_t>(grid_map.width),
                                 static_cast<size_t>(grid_map.height),
                                 grid_map.data.size(),
                                 limits,
                                 estimate,
                                 result.error)) {
    return result;
  }
  (void)estimate;
  return finishCreation(
    Polymap(grid_map, start_x, start_y, goal_x, goal_y, start_position, goal_position, stop_token));
}

PolymapCreateResult Polymap::create(const GridMap& grid_map,
                                    int start_x,
                                    int start_y,
                                    const Point2d& start_position,
                                    const std::vector<PolymapEndpoint>& goals,
                                    const StopToken& stop_token,
                                    const PlanningLimits& limits) {
  PolymapCreateResult result;
  MapResourceEstimate estimate;
  if (!validateMapResourceBudget(static_cast<size_t>(grid_map.width),
                                 static_cast<size_t>(grid_map.height),
                                 grid_map.data.size(),
                                 limits,
                                 estimate,
                                 result.error)) {
    return result;
  }
  (void)estimate;
  return finishCreation(Polymap(grid_map, start_x, start_y, start_position, goals, stop_token));
}

PolymapCreateResult Polymap::createForReferenceShortening(
  const GridMap& grid_map,
  int start_x,
  int start_y,
  const Point2d& start_position,
  const std::vector<PolymapEndpoint>& goals,
  const StopToken& stop_token,
  const PlanningLimits& limits) {
  PolymapCreateResult result;
  MapResourceEstimate estimate;
  if (!validateMapResourceBudget(static_cast<size_t>(grid_map.width),
                                 static_cast<size_t>(grid_map.height),
                                 grid_map.data.size(),
                                 limits,
                                 estimate,
                                 result.error)) {
    return result;
  }
  const size_t remaining_map_bytes = limits.max_map_bytes - estimate.estimated_bytes;
  const size_t max_raw_contour_vertices =
    remaining_map_bytes / kEstimatedReferenceTopologyBytesPerRawContourVertex;
  return finishCreation(
    Polymap(grid_map,
            start_x,
            start_y,
            start_position,
            goals,
            stop_token,
            false,
            RawContourResourceBudget{max_raw_contour_vertices, limits.max_map_bytes}));
}

Polymap::Polymap(const GridMap& grid_map,
                 int start_x,
                 int start_y,
                 int goal_x,
                 int goal_y,
                 const Point2d& start_position,
                 const Point2d& goal_position,
                 const StopToken& stop_token)
  : Polymap(grid_map,
            start_x,
            start_y,
            start_position,
            std::vector<PolymapEndpoint>{{goal_x, goal_y, goal_position}},
            stop_token) {}

Polymap::Polymap(const GridMap& grid_map,
                 int start_x,
                 int start_y,
                 const Point2d& start_position,
                 const std::vector<PolymapEndpoint>& goals,
                 const StopToken& stop_token)
  : Polymap(grid_map, start_x, start_y, start_position, goals, stop_token, true, std::nullopt) {}

Polymap::Polymap(const GridMap& grid_map,
                 int start_x,
                 int start_y,
                 const Point2d& start_position,
                 const std::vector<PolymapEndpoint>& goals,
                 const StopToken& stop_token,
                 bool simplify_obstacle_contours,
                 std::optional<RawContourResourceBudget> raw_contour_budget)
  : xsize_(0), ysize_(0) {
  const auto stop_construction = [&]() { clearStoppedConstructionState(); };
  const auto reject = [&](const std::string& message) { construction_error_ = message; };

  if (stop_token.poll()) {
    stop_construction();
    return;
  }

  size_t cell_count = 0;
  if (!prepareGridStorage(grid_map, xsize_, ysize_, cell_count, construction_error_)) {
    return;
  }
  data_ = grid_map.data;
  vertices_location_x_flat_.resize(cell_count, -1);
  vertices_location_y_flat_.resize(cell_count, -1);

  if (!std::isfinite(start_position.first) || !std::isfinite(start_position.second)) {
    reject("Start position is not finite");
    return;
  }
  if (goals.empty()) {
    reject("At least one goal is required");
    return;
  }
  if (std::floor(start_position.first) != static_cast<double>(start_x) ||
      std::floor(start_position.second) != static_cast<double>(start_y)) {
    reject("Start continuous position does not belong to the supplied grid cell");
    return;
  }
  if (start_x < 0 || start_y < 0 || start_x >= xsize_ || start_y >= ysize_) {
    reject("Start grid cell is outside map bounds");
    return;
  }
  const size_t start_index =
    static_cast<size_t>(start_x) + static_cast<size_t>(start_y) * static_cast<size_t>(xsize_);
  if (start_index >= data_.size() || data_[start_index] != 0) {
    reject("Start grid cell is occupied");
    return;
  }

  for (const auto& goal : goals) {
    if (!std::isfinite(goal.position.first) || !std::isfinite(goal.position.second)) {
      reject("Goal position is not finite");
      return;
    }
    if (std::floor(goal.position.first) != static_cast<double>(goal.cell_x) ||
        std::floor(goal.position.second) != static_cast<double>(goal.cell_y)) {
      reject("Goal continuous position does not belong to the supplied grid cell");
      return;
    }
    if (goal.cell_x < 0 || goal.cell_y < 0 || goal.cell_x >= xsize_ || goal.cell_y >= ysize_) {
      reject("Goal grid cell is outside map bounds");
      return;
    }
    const size_t goal_index = static_cast<size_t>(goal.cell_x) +
                              static_cast<size_t>(goal.cell_y) * static_cast<size_t>(xsize_);
    if (goal_index >= data_.size() || data_[goal_index] != 0) {
      reject("Goal grid cell is occupied");
      return;
    }
  }

  const std::optional<size_t> max_raw_contour_vertices =
    raw_contour_budget ? std::optional<size_t>(raw_contour_budget->max_vertices) : std::nullopt;
  const OperationStatus obstacle_status =
    getPolyObstacles(start_x, start_y, goals, stop_token, max_raw_contour_vertices);
  if (obstacle_status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  solution_exist_ = obstacle_status == OperationStatus::success;
  if (!solution_exist_) {
    if (!construction_error_.empty()) {
      if (raw_contour_budget) {
        construction_error_ +=
          " derived from max_map_bytes=" +
          std::to_string(raw_contour_budget->max_map_bytes);
      }
      return;
    }
    no_path_ = true;
    construction_error_ = "Start and every goal must be in the same reachable free-space component";
    return;
  }

  if (!simplify_obstacle_contours && raw_contour_budget) {
    const size_t max_raw_contour_vertices = raw_contour_budget->max_vertices;
    size_t raw_contour_vertices = 0;
    for (const auto& obstacle : obs_) {
      if (stop_token.poll()) {
        stop_construction();
        return;
      }
      if (raw_contour_vertices > max_raw_contour_vertices ||
          obstacle.ordered_vertices_.size() >
            max_raw_contour_vertices - raw_contour_vertices) {
        construction_error_ =
          "Reference-shortening raw contour exceeds the vertex budget " +
          std::to_string(max_raw_contour_vertices) + " derived from max_map_bytes=" +
          std::to_string(raw_contour_budget->max_map_bytes);
        return;
      }
      raw_contour_vertices += obstacle.ordered_vertices_.size();
    }
  }

  OperationStatus status = OperationStatus::success;
  const bool has_outer_contour = hasClockwiseOuterContour(obs_, stop_token);
  if (stop_token.poll()) {
    stop_construction();
    return;
  }

  // Reject continuous endpoints against the unsimplified contours first.  In
  // particular, a point on a raw contour must never become accepted merely
  // because a later chord removes the incident contour vertices.
  std::string endpoint_error;
  if (has_outer_contour) {
    status = validateFreeSpaceInterior(start_position, stop_token, &endpoint_error);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Invalid start position: " + endpoint_error;
      return;
    }
    for (const auto& goal : goals) {
      status = validateFreeSpaceInterior(goal.position, stop_token, &endpoint_error);
      if (status == OperationStatus::stopped) {
        stop_construction();
        return;
      }
      if (status == OperationStatus::failure) {
        construction_error_ = "Invalid goal position: " + endpoint_error;
        return;
      }
    }
  }

  status = validateObstacleTopology(construction_error_, false, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;

  if (simplify_obstacle_contours) {
    std::vector<Point2d> protected_points;
    protected_points.reserve(goals.size() + 1);
    protected_points.emplace_back(start_position);
    for (const auto& goal : goals) protected_points.emplace_back(goal.position);
    status = simplifyPolyObstaclesImpl(protected_points, stop_token)
               ? OperationStatus::success
               : (stop_token.poll() ? OperationStatus::stopped : OperationStatus::failure);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure) {
      construction_error_ = "Obstacle simplification failed";
      return;
    }

    status = validateObstacleTopology(construction_error_, true, stop_token);
    if (status == OperationStatus::stopped) {
      stop_construction();
      return;
    }
    if (status == OperationStatus::failure)
      return;

    if (has_outer_contour) {
      status = validateFreeSpaceInterior(start_position, stop_token, &endpoint_error);
      if (status == OperationStatus::stopped) {
        stop_construction();
        return;
      }
      if (status == OperationStatus::failure) {
        construction_error_ = "Invalid start position: " + endpoint_error;
        return;
      }
      for (const auto& goal : goals) {
        status = validateFreeSpaceInterior(goal.position, stop_token, &endpoint_error);
        if (status == OperationStatus::stopped) {
          stop_construction();
          return;
        }
        if (status == OperationStatus::failure) {
          construction_error_ = "Invalid goal position: " + endpoint_error;
          return;
        }
      }
    }
  }

  status = registerVertices(construction_error_, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;

  status = constructCGALRelated(&validateCDTWithCGAL, construction_error_, stop_token);
  if (status == OperationStatus::stopped) {
    stop_construction();
    return;
  }
  if (status == OperationStatus::failure)
    return;
}

bool Polymap::registerVertices(std::string& error) {
  const StopToken no_stop;
  return registerVerticesImpl(error, no_stop);
}

OperationStatus Polymap::registerVertices(std::string& error, const StopToken& stop_token) {
  error.clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (registerVerticesImpl(error, stop_token))
    return OperationStatus::success;
  if (stop_token.poll()) {
    error.clear();
    return OperationStatus::stopped;
  }
  return OperationStatus::failure;
}

bool Polymap::registerVerticesImpl(std::string& error, const StopToken& stop_token) {
  error.clear();
  if (stop_token.poll())
    return false;
  const size_t registry_size = static_cast<size_t>(xsize_) * static_cast<size_t>(ysize_);
  std::vector<int> candidate_obstacles(registry_size, -1);
  std::vector<int> candidate_vertices(registry_size, -1);

  for (size_t i = 0; i < obs_.size(); ++i) {
    if (stop_token.poll())
      return false;
    for (size_t j = 0; j < obs_[i].ordered_vertices_.size(); ++j) {
      if (stop_token.poll())
        return false;
      const int x = obs_[i].ordered_vertices_[j].first;
      const int y = obs_[i].ordered_vertices_[j].second;
      if (x < 0 || y < 0 || x >= xsize_ || y >= ysize_) {
        error = "Unsupported obstacle topology: obstacle " + std::to_string(i) + " vertex " +
                std::to_string(j) + " at (" + std::to_string(x) + ", " + std::to_string(y) +
                ") is outside the vertex registry";
        return false;
      }

      const size_t index =
        static_cast<size_t>(y) * static_cast<size_t>(xsize_) + static_cast<size_t>(x);
      if (candidate_obstacles[index] >= 0) {
        error = "Unsupported obstacle topology: obstacle " + std::to_string(i) + " vertex " +
                std::to_string(j) + " would overwrite obstacle " +
                std::to_string(candidate_obstacles[index]) + " vertex " +
                std::to_string(candidate_vertices[index]) + " at shared vertex (" +
                std::to_string(x) + ", " + std::to_string(y) + ")";
        return false;
      }
      candidate_obstacles[index] = static_cast<int>(i);
      candidate_vertices[index] = static_cast<int>(j);
    }
  }

  if (stop_token.poll())
    return false;
  vertices_location_x_flat_.swap(candidate_obstacles);
  vertices_location_y_flat_.swap(candidate_vertices);
  return true;
}

void Polymap::simplifyPolyObstacles(int start_x, int start_y, int goal_x, int goal_y) {
  const StopToken no_stop;
  (void)simplifyPolyObstaclesImpl(
    Point2d{static_cast<double>(start_x), static_cast<double>(start_y)},
    Point2d{static_cast<double>(goal_x), static_cast<double>(goal_y)},
    no_stop);
}

OperationStatus Polymap::simplifyPolyObstacles(
  int start_x, int start_y, int goal_x, int goal_y, const StopToken& stop_token) {
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (simplifyPolyObstaclesImpl(Point2d{static_cast<double>(start_x), static_cast<double>(start_y)},
                                Point2d{static_cast<double>(goal_x), static_cast<double>(goal_y)},
                                stop_token)) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    return OperationStatus::success;
  }
  return stop_token.poll() ? OperationStatus::stopped : OperationStatus::failure;
}

OperationStatus Polymap::simplifyPolyObstacles(const Point2d& start,
                                               const Point2d& goal,
                                               const StopToken& stop_token) {
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (simplifyPolyObstaclesImpl(start, goal, stop_token)) {
    if (stop_token.poll())
      return OperationStatus::stopped;
    return OperationStatus::success;
  }
  return stop_token.poll() ? OperationStatus::stopped : OperationStatus::failure;
}

void Polymap::simplifyPolyObstacles(const Point2d& start, const Point2d& goal) {
  const StopToken no_stop;
  (void)simplifyPolyObstaclesImpl(start, goal, no_stop);
}

bool Polymap::simplifyPolyObstaclesImpl(const Point2d& start,
                                        const Point2d& goal,
                                        const StopToken& stop_token) {
  return simplifyPolyObstaclesImpl(std::vector<Point2d>{start, goal}, stop_token);
}

bool Polymap::simplifyPolyObstaclesImpl(const std::vector<Point2d>& protected_points,
                                        const StopToken& stop_token) {
  for (auto iter = obs_.begin(); iter != obs_.end(); ++iter) {
    if (stop_token.poll())
      return false;
    int prev, curr, next;
    bool stable = false;
    curr = 0;
    bool simplifable;
    int x1, y1, x2, y2, x3, y3;

    while (1) {
      if (stop_token.poll())
        return false;
      if (iter->ordered_vertices_.size() <= 3)
        break;
      prev = (curr - 1 + static_cast<int>(iter->ordered_vertices_.size())) %
             static_cast<int>(iter->ordered_vertices_.size());
      next = (curr + 1) % iter->ordered_vertices_.size();

      x1 = iter->ordered_vertices_[prev].first;
      y1 = iter->ordered_vertices_[prev].second;
      x2 = iter->ordered_vertices_[curr].first;
      y2 = iter->ordered_vertices_[curr].second;
      x3 = iter->ordered_vertices_[next].first;
      y3 = iter->ordered_vertices_[next].second;

      const exact_geometry::Point previous_point(x1, y1);
      const exact_geometry::Point current_point(x2, y2);
      const exact_geometry::Point next_point(x3, y3);
      const CGAL::Orientation turn = CGAL::orientation(previous_point, current_point, next_point);

      if (turn == CGAL::COLLINEAR) {
        // Only a true middle point is redundant. A collinear backtracking
        // spike changes the contour and must not be erased.
        simplifable =
          exact_geometry::isRemovableCollinearMiddle(previous_point, current_point, next_point);
      } else {
        simplifable = true;
        if (turn == CGAL::RIGHT_TURN) {
          for (const auto& point : protected_points) {
            if (isInTri(x1, y1, x2, y2, x3, y3, point.first, point.second)) {
              simplifable = false;
              break;
            }
          }

          double testx, testy;
          if (simplifable) {
            for (auto iter2 = obs_.begin(); iter2 != obs_.end(); ++iter2) {
              if (stop_token.poll())
                return false;
              for (auto iter3 = iter2->ordered_vertices_.begin();
                   iter3 != iter2->ordered_vertices_.end();
                   ++iter3) {
                if (stop_token.poll())
                  return false;
                testx = iter3->first;
                testy = iter3->second;
                if (isInTri(x1, y1, x2, y2, x3, y3, testx, testy)) {
                  if (iter2 - obs_.begin() != iter - obs_.begin())
                    simplifable = false;
                  else {
                    if (iter3 - iter2->ordered_vertices_.begin() != prev &&
                        iter3 - iter2->ordered_vertices_.begin() != curr &&
                        iter3 - iter2->ordered_vertices_.begin() != next)
                      simplifable = false;
                  }
                }
              }
            }
          }
        } else {
          simplifable = false;
        }
      }

      if (simplifable) {
        const size_t obstacle_index = static_cast<size_t>(iter - obs_.begin());
        simplifable = simplificationChordIsTopologicallySafeImpl(
          obstacle_index, static_cast<size_t>(curr), stop_token);
        if (stop_token.poll())
          return false;
      }

      if (simplifable) {
        iter->ordered_vertices_.erase(iter->ordered_vertices_.begin() + curr);
        if (curr >= static_cast<int>(iter->ordered_vertices_.size()))
          curr = static_cast<int>(iter->ordered_vertices_.size()) - 1;
        stable = false;
        continue;
      } else {
        if (curr == 0) {
          if (!stable)
            stable = true;
          else
            break;
        }
      }

      curr++;
      if (curr >= static_cast<int>(iter->ordered_vertices_.size()))
        curr = 0;
    }
  }
  return true;
}

}  // namespace raystar
