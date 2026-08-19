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

bool validateReducedDirectedPortalWitness(const TriangleCorridor& corridor, std::string* error) {
  if (error)
    error->clear();
  const auto fail = [&](const std::string& message) {
    if (error)
      *error = message;
    return false;
  };

  if (corridor.triangle_occurrences.empty())
    return fail("A lifted-sleeve witness requires a containing triangle occurrence");
  if (corridor.portals.size() == std::numeric_limits<size_t>::max() ||
      corridor.triangle_occurrences.size() != corridor.portals.size() + 1) {
    return fail("Triangle-occurrence and directed-portal cardinalities are inconsistent");
  }

  for (size_t index = 0; index < corridor.portals.size(); ++index) {
    const auto& portal = corridor.portals[index];
    if (portal.from_triangle != corridor.triangle_occurrences[index] ||
        portal.to_triangle != corridor.triangle_occurrences[index + 1]) {
      return fail("Directed portal " + std::to_string(index) +
                  " does not bind its adjacent triangle occurrences");
    }
    if (portal.from_triangle == portal.to_triangle) {
      return fail("Directed portal " + std::to_string(index) + " is a triangle self-transition");
    }
    if (!std::isfinite(portal.left.first) || !std::isfinite(portal.left.second) ||
        !std::isfinite(portal.right.first) || !std::isfinite(portal.right.second)) {
      return fail("Directed portal " + std::to_string(index) + " contains a non-finite endpoint");
    }
    if (portal.left == portal.right) {
      return fail("Directed portal " + std::to_string(index) + " has zero geometric width");
    }
  }

  // traceFreeSpacePath() reduces the dual walk as a stack.  Retain later
  // repeated occurrences (which encode winding), but reject an adjacent
  // portal followed immediately by its inverse A->B->A.
  for (size_t index = 2; index < corridor.triangle_occurrences.size(); ++index) {
    if (corridor.triangle_occurrences[index] == corridor.triangle_occurrences[index - 2]) {
      return fail("Lifted-sleeve witness contains an unreduced immediate portal reversal at " +
                  std::to_string(index - 2));
    }
  }
  return true;
}

bool sameReducedDirectedPortalWitness(const TriangleCorridor& reference,
                                      const TriangleCorridor& candidate,
                                      std::string* error) {
  if (error)
    error->clear();
  const auto fail = [&](const std::string& message) {
    if (error)
      *error = message;
    return false;
  };

  std::string validation_error;
  if (!validateReducedDirectedPortalWitness(reference, &validation_error))
    return fail("Reference lifted-sleeve witness is invalid: " + validation_error);
  if (!validateReducedDirectedPortalWitness(candidate, &validation_error))
    return fail("Candidate lifted-sleeve witness is invalid: " + validation_error);

  if (reference.triangle_occurrences.size() != candidate.triangle_occurrences.size()) {
    return fail("Lifted-sleeve triangle-occurrence counts differ");
  }
  for (size_t index = 0; index < reference.triangle_occurrences.size(); ++index) {
    if (reference.triangle_occurrences[index] != candidate.triangle_occurrences[index]) {
      return fail("Lifted-sleeve triangle occurrence differs at ordinal " + std::to_string(index));
    }
  }

  // Cardinalities were validated above and the equal triangle counts imply
  // equal portal counts.  Compare each occurrence in order; never sort or
  // deduplicate repeated portal geometry.
  for (size_t index = 0; index < reference.portals.size(); ++index) {
    const auto& expected = reference.portals[index];
    const auto& actual = candidate.portals[index];
    if (expected.from_triangle != actual.from_triangle ||
        expected.to_triangle != actual.to_triangle) {
      return fail("Directed portal identity differs at occurrence " + std::to_string(index));
    }
    if (expected.left != actual.left || expected.right != actual.right) {
      return fail("Directed portal geometry differs at occurrence " + std::to_string(index));
    }
  }
  return true;
}

void Polymap::clearCGALRelatedState() {
  cdt_ready_ = false;
  cdt_.clear();
  facets_.clear();
  triangle_faces_.clear();
  triangle_edges_.clear();
  cdt_table_.clear();
  cdt_ver_num_ = 0;
  visibility_storage_.clear();
}

void Polymap::clearStoppedConstructionState() {
  solution_exist_ = false;
  no_path_ = false;
  construction_stopped_ = true;
  construction_error_.clear();
  obs_.clear();
  std::fill(vertices_location_x_flat_.begin(), vertices_location_x_flat_.end(), -1);
  std::fill(vertices_location_y_flat_.begin(), vertices_location_y_flat_.end(), -1);
  clearCGALRelatedState();
}

bool Polymap::constructCGALRelated(CdtValidator validator, std::string& error) {
  const StopToken no_stop;
  return constructCGALRelatedImpl(validator, error, no_stop);
}

OperationStatus Polymap::constructCGALRelated(CdtValidator validator,
                                              std::string& error,
                                              const StopToken& stop_token) {
  error.clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  if (constructCGALRelatedImpl(validator, error, stop_token))
    return OperationStatus::success;
  if (stop_token.poll()) {
    error.clear();
    return OperationStatus::stopped;
  }
  return OperationStatus::failure;
}

bool Polymap::constructCGALRelatedImpl(CdtValidator validator,
                                       std::string& error,
                                       const StopToken& stop_token) {
  using CDT = constrained_delaunay_triangulation::CDT;

  error.clear();
  clearCGALRelatedState();

  const auto fail = [&](std::string message) {
    clearCGALRelatedState();
    error = std::move(message);
    return false;
  };

  if (validator == nullptr)
    return fail("CDT validation function is unavailable");
  if (stop_token.poll())
    return false;

  try {
    CDT candidate_cdt;
    for (const auto& ob : obs_) {
      if (stop_token.poll())
        return false;
      for (auto it = ob.ordered_vertices_.begin(); it != ob.ordered_vertices_.end(); ++it) {
        if (stop_token.poll())
          return false;
        auto next = std::next(it);
        if (next == ob.ordered_vertices_.end())
          next = ob.ordered_vertices_.begin();
        candidate_cdt.insert_constraint(
          constrained_delaunay_triangulation::Point(it->first, it->second),
          constrained_delaunay_triangulation::Point(next->first, next->second));
        if (stop_token.poll())
          return false;
      }
    }

    if (stop_token.poll())
      return false;
    const bool valid_cdt = validator(candidate_cdt);
    if (stop_token.poll())
      return false;
    if (!valid_cdt)
      return fail("CGAL reported that the constrained triangulation is invalid");

    std::vector<std::vector<std::pair<int, int>>> candidate_facets;
    std::unordered_map<long long, int> candidate_table;
    const int candidate_vertex_count = static_cast<int>(candidate_cdt.number_of_vertices());
    const long long map_area = static_cast<long long>(xsize_) * static_cast<long long>(ysize_);
    const auto vertex_key = [&](const std::pair<int, int>& vertex) {
      return static_cast<long long>(vertex.first) +
             static_cast<long long>(vertex.second) * static_cast<long long>(xsize_);
    };
    const auto directed_edge_key = [&](const std::pair<int, int>& from,
                                       const std::pair<int, int>& to) {
      return vertex_key(from) + vertex_key(to) * map_area;
    };

    int count = 0;
    for (auto fit = candidate_cdt.finite_faces_begin(); fit != candidate_cdt.finite_faces_end();
         ++fit) {
      if (stop_token.poll())
        return false;
      candidate_facets.emplace_back(
        std::vector<std::pair<int, int>>{{static_cast<int>(fit->vertex(0)->point().x()),
                                          static_cast<int>(fit->vertex(0)->point().y())},
                                         {static_cast<int>(fit->vertex(1)->point().x()),
                                          static_cast<int>(fit->vertex(1)->point().y())},
                                         {static_cast<int>(fit->vertex(2)->point().x()),
                                          static_cast<int>(fit->vertex(2)->point().y())}});
      const auto& facet = candidate_facets.back();
      fit->info().stable_id = count;
      fit->info().is_free = false;
      candidate_table[directed_edge_key(facet[0], facet[1])] = count;
      candidate_table[directed_edge_key(facet[1], facet[2])] = count;
      candidate_table[directed_edge_key(facet[2], facet[0])] = count;
      ++count;
    }

    if (stop_token.poll())
      return false;
    cdt_ = std::move(candidate_cdt);
    facets_ = std::move(candidate_facets);
    cdt_table_ = std::move(candidate_table);
    cdt_ver_num_ = candidate_vertex_count;
    cdt_ready_ = true;
    const OperationStatus triangle_environment_status =
      buildTriangleEnvironment(error, stop_token);
    if (triangle_environment_status == OperationStatus::stopped)
      return false;
    if (triangle_environment_status == OperationStatus::failure)
      return fail(error.empty() ? "Could not build the free-triangle environment" : error);
    return true;
  } catch (const CDT::Intersection_of_constraints_exception&) {
    return fail("Obstacle constraints intersect, overlap, or form an unsupported T-junction");
  } catch (const CGAL::Failure_exception& exception) {
    return fail("CGAL rejected the constrained triangulation: " + std::string(exception.what()));
  } catch (const std::exception& exception) {
    return fail("Unexpected error while constructing the CDT: " + std::string(exception.what()));
  } catch (...) {
    return fail("Unknown error while constructing the CDT");
  }
}

OperationStatus Polymap::buildTriangleEnvironment(std::string& error,
                                                  const StopToken& stop_token) {
  using EdgeKey = std::pair<Point2d, Point2d>;
  struct EdgeUse {
    int face = -1;
    int local_edge = -1;
  };
  struct EdgeAccumulator {
    std::vector<EdgeUse> uses;
    bool constrained = false;
  };

  error.clear();
  triangle_faces_.clear();
  triangle_edges_.clear();
  const auto stopped = [&]() {
    error.clear();
    triangle_faces_.clear();
    triangle_edges_.clear();
    return OperationStatus::stopped;
  };
  const auto fail = [&](std::string message) {
    error = std::move(message);
    triangle_faces_.clear();
    triangle_edges_.clear();
    return OperationStatus::failure;
  };
  if (stop_token.poll())
    return stopped();
  if (!cdt_ready_ || facets_.empty() || obs_.empty()) {
    return fail("Triangle environment requires a ready CDT and a reachable outer contour");
  }
  if (facets_.size() > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
    return fail("Triangle count exceeds the stable TriangleId range");
  }

  const auto canonical_edge = [](Point2d first, Point2d second) {
    if (second < first)
      std::swap(first, second);
    return EdgeKey{first, second};
  };

  std::set<EdgeKey> constrained_edges;
  for (const auto& obstacle : obs_) {
    if (stop_token.poll())
      return stopped();
    const auto& vertices = obstacle.ordered_vertices_;
    for (size_t index = 0; index < vertices.size(); ++index) {
      if (stop_token.poll())
        return stopped();
      const auto& first = vertices[index];
      const auto& second = vertices[(index + 1) % vertices.size()];
      constrained_edges.insert(
        canonical_edge({static_cast<double>(first.first), static_cast<double>(first.second)},
                       {static_cast<double>(second.first), static_cast<double>(second.second)}));
    }
  }

  size_t outer_index = 0;
  long double outer_area = 0.0L;
  bool found_clockwise_outer = false;
  for (size_t obstacle_index = 0; obstacle_index < obs_.size(); ++obstacle_index) {
    if (stop_token.poll())
      return stopped();
    const auto& vertices = obs_[obstacle_index].ordered_vertices_;
    long double twice_area = 0.0L;
    for (size_t index = 0; index < vertices.size(); ++index) {
      if (stop_token.poll())
        return stopped();
      const auto& from = vertices[index];
      const auto& to = vertices[(index + 1) % vertices.size()];
      twice_area += static_cast<long double>(from.first) * to.second -
                    static_cast<long double>(to.first) * from.second;
    }
    if ((!found_clockwise_outer && twice_area < 0.0L) ||
        (found_clockwise_outer && twice_area < outer_area)) {
      outer_index = obstacle_index;
      outer_area = twice_area;
      found_clockwise_outer = true;
    } else if (!found_clockwise_outer && (obstacle_index == 0 || twice_area > outer_area)) {
      outer_index = obstacle_index;
      outer_area = twice_area;
    }
  }

  triangle_faces_.resize(facets_.size());
  std::map<EdgeKey, EdgeAccumulator> edges;
  for (size_t face_index = 0; face_index < facets_.size(); ++face_index) {
    if (stop_token.poll())
      return stopped();
    const auto& facet = facets_[face_index];
    if (facet.size() != 3) {
      return fail("CDT facet does not contain exactly three vertices");
    }
    auto& face = triangle_faces_[face_index];
    long double centroid_x = 0.0L;
    long double centroid_y = 0.0L;
    for (size_t vertex_index = 0; vertex_index < 3; ++vertex_index) {
      face.vertices[vertex_index] = {static_cast<double>(facet[vertex_index].first),
                                     static_cast<double>(facet[vertex_index].second)};
      centroid_x += facet[vertex_index].first;
      centroid_y += facet[vertex_index].second;
    }
    centroid_x /= 3.0L;
    centroid_y /= 3.0L;
    bool point_is_inside = false;
    if (pointInPolygon(obs_[outer_index].ordered_vertices_,
                       centroid_x,
                       centroid_y,
                       stop_token,
                       point_is_inside) == OperationStatus::stopped) {
      return stopped();
    }
    face.is_free = point_is_inside;
    for (size_t obstacle_index = 0; obstacle_index < obs_.size() && face.is_free;
         ++obstacle_index) {
      if (stop_token.poll())
        return stopped();
      if (obstacle_index == outer_index)
        continue;
      if (pointInPolygon(obs_[obstacle_index].ordered_vertices_,
                         centroid_x,
                         centroid_y,
                         stop_token,
                         point_is_inside) == OperationStatus::stopped) {
        return stopped();
      }
      if (point_is_inside)
        face.is_free = false;
    }

    for (size_t edge_index = 0; edge_index < 3; ++edge_index) {
      if (stop_token.poll())
        return stopped();
      const EdgeKey key =
        canonical_edge(face.vertices[edge_index], face.vertices[(edge_index + 1) % 3]);
      auto& accumulator = edges[key];
      accumulator.uses.push_back({static_cast<int>(face_index), static_cast<int>(edge_index)});
      accumulator.constrained = constrained_edges.find(key) != constrained_edges.end();
    }
  }

  for (auto& [key, accumulator] : edges) {
    if (stop_token.poll())
      return stopped();
    if (accumulator.uses.empty() || accumulator.uses.size() > 2) {
      return fail("CDT edge has an invalid number of incident finite faces");
    }
    TriangleMeshEdge edge;
    edge.a = key.first;
    edge.b = key.second;
    edge.constrained = accumulator.constrained;
    for (size_t use_index = 0; use_index < accumulator.uses.size(); ++use_index) {
      const auto& use = accumulator.uses[use_index];
      edge.faces[use_index] = use.face;
      triangle_faces_[static_cast<size_t>(use.face)]
        .constrained[static_cast<size_t>(use.local_edge)] = edge.constrained;
    }
    if (accumulator.uses.size() == 2 && !edge.constrained) {
      const auto& first = accumulator.uses[0];
      const auto& second = accumulator.uses[1];
      if (triangle_faces_[static_cast<size_t>(first.face)].is_free &&
          triangle_faces_[static_cast<size_t>(second.face)].is_free) {
        triangle_faces_[static_cast<size_t>(first.face)]
          .neighbors[static_cast<size_t>(first.local_edge)] = second.face;
        triangle_faces_[static_cast<size_t>(second.face)]
          .neighbors[static_cast<size_t>(second.local_edge)] = first.face;
      }
    }
    triangle_edges_.push_back(edge);
  }

  size_t free_count = 0;
  for (auto fit = cdt_.finite_faces_begin(); fit != cdt_.finite_faces_end(); ++fit) {
    if (stop_token.poll())
      return stopped();
    const int id = fit->info().stable_id;
    if (id < 0 || id >= static_cast<int>(triangle_faces_.size())) {
      return fail("CDT face lost its stable ID while building the triangle environment");
    }
    fit->info().is_free = triangle_faces_[static_cast<size_t>(id)].is_free;
    if (fit->info().is_free)
      ++free_count;
  }
  if (free_count == 0) {
    return fail("Triangle environment contains no free faces");
  }
  if (stop_token.poll())
    return stopped();
  return OperationStatus::success;
}

size_t Polymap::freeTriangleCount() const noexcept {
  return static_cast<size_t>(std::count_if(
    triangle_faces_.begin(), triangle_faces_.end(), [](const auto& face) { return face.is_free; }));
}

OperationStatus Polymap::traceFreeSpacePath(const std::vector<Point2d>& path,
                                            TriangleCorridor& corridor,
                                            const StopToken& stop_token,
                                            std::string* error) const {
  corridor = {};
  if (error)
    error->clear();
  if (stop_token.poll())
    return OperationStatus::stopped;
  return traceFreeSpacePathImpl(path, corridor, stop_token, error);
}

OperationStatus Polymap::traceFreeSpacePathImpl(const std::vector<Point2d>& input_path,
                                                TriangleCorridor& corridor,
                                                const StopToken& stop_token,
                                                std::string* error) const {
  const auto stopped = [&]() {
    corridor = {};
    if (error)
      error->clear();
    return OperationStatus::stopped;
  };
  const auto fail = [&](const std::string& message) {
    corridor = {};
    if (error)
      *error = message;
    return OperationStatus::failure;
  };
  if (!cdt_ready_ || triangle_faces_.empty())
    return fail("Free-triangle environment is not ready");
  if (input_path.empty())
    return fail("A reference polyline requires at least one point");

  std::vector<Point2d> path;
  path.reserve(input_path.size());
  for (const auto& point : input_path) {
    if (!std::isfinite(point.first) || !std::isfinite(point.second))
      return fail("Reference polyline contains a non-finite point");
    if (path.empty() || !samePoint(path.back(), point))
      path.push_back(point);
  }
  const auto locate_candidates = [&](const Point2d& point) {
    std::vector<int> candidates;
    using CDT = constrained_delaunay_triangulation::CDT;
    CDT::Locate_type locate_type;
    int local_index = -1;
    const auto face =
      cdt_.locate(constrained_delaunay_triangulation::Point(point.first, point.second),
                  locate_type,
                  local_index);
    const auto add_face = [&](const CDT::Face_handle& candidate) {
      if (candidate == CDT::Face_handle() || cdt_.is_infinite(candidate))
        return;
      const int id = candidate->info().stable_id;
      if (id >= 0 && id < static_cast<int>(triangle_faces_.size()) &&
          triangle_faces_[static_cast<size_t>(id)].is_free) {
        candidates.push_back(id);
      }
    };
    if (locate_type == CDT::FACE) {
      add_face(face);
    } else if (locate_type == CDT::EDGE) {
      add_face(face);
      add_face(face->neighbor(local_index));
    } else if (locate_type == CDT::VERTEX) {
      const auto vertex = face->vertex(local_index);
      auto incident = cdt_.incident_faces(vertex, face);
      const auto done = incident;
      if (incident != 0) {
        do {
          add_face(incident);
          ++incident;
        } while (incident != done);
      }
    }
    std::sort(candidates.begin(), candidates.end());
    candidates.erase(std::unique(candidates.begin(), candidates.end()), candidates.end());
    return candidates;
  };

  // A constant path is the identity transition. Keep one containing triangle
  // as its zero-portal sleeve so callers can distinguish it from an invalid
  // or out-of-free-space reference.
  if (path.size() == 1) {
    const auto candidates = locate_candidates(path.front());
    if (candidates.empty())
      return fail("The constant reference point is outside the free-triangle environment");
    corridor.triangle_occurrences.push_back(static_cast<uint32_t>(candidates.front()));
    return OperationStatus::success;
  }

  const auto connect_at_point = [&](int from, int to, const Point2d& point) {
    std::vector<int> empty;
    if (from == to)
      return std::vector<int>{from};
    std::vector<int> parent(triangle_faces_.size(), -2);
    std::queue<int> pending;
    parent[static_cast<size_t>(from)] = -1;
    pending.push(from);
    while (!pending.empty()) {
      if (stop_token.poll())
        return empty;
      const int current = pending.front();
      pending.pop();
      std::vector<int> neighbors;
      for (const int neighbor : triangle_faces_[static_cast<size_t>(current)].neighbors) {
        if (neighbor >= 0)
          neighbors.push_back(neighbor);
      }
      std::sort(neighbors.begin(), neighbors.end());
      neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
      for (const int neighbor : neighbors) {
        if (parent[static_cast<size_t>(neighbor)] != -2 ||
            !pointInTriangleClosed(triangle_faces_[static_cast<size_t>(neighbor)].vertices,
                                   point)) {
          continue;
        }
        parent[static_cast<size_t>(neighbor)] = current;
        if (neighbor == to) {
          std::vector<int> connection;
          for (int trace = to; trace >= 0; trace = parent[static_cast<size_t>(trace)])
            connection.push_back(trace);
          std::reverse(connection.begin(), connection.end());
          return connection;
        }
        pending.push(neighbor);
      }
    }
    return empty;
  };

  std::string segment_failure;
  const auto segment_faces = [&](const Point2d& start, const Point2d& goal) {
    segment_failure.clear();
    std::vector<int> faces;
    std::vector<long double> events{0.0L, 1.0L};
    events.reserve(triangle_edges_.size() + 2);
    for (const auto& edge : triangle_edges_) {
      if (stop_token.poll())
        return faces;
      appendSegmentEdgeEvents(start, goal, edge.a, edge.b, events);
    }
    std::sort(events.begin(), events.end());
    constexpr long double event_tolerance = 1.0e-13L;
    events.erase(std::unique(events.begin(),
                             events.end(),
                             [](long double first, long double second) {
                               return std::abs(first - second) <= event_tolerance;
                             }),
                 events.end());

    for (size_t interval = 0; interval + 1 < events.size(); ++interval) {
      if (stop_token.poll())
        return std::vector<int>{};
      if (events[interval + 1] - events[interval] <= event_tolerance)
        continue;
      const long double parameter = (events[interval] + events[interval + 1]) / 2.0L;
      const Point2d midpoint{
        static_cast<double>(static_cast<long double>(start.first) +
                            parameter * (static_cast<long double>(goal.first) - start.first)),
        static_cast<double>(static_cast<long double>(start.second) +
                            parameter * (static_cast<long double>(goal.second) - start.second))};
      const auto candidates = locate_candidates(midpoint);
      if (candidates.empty()) {
        std::ostringstream stream;
        stream << "no free face contains interval midpoint (" << midpoint.first << ", "
               << midpoint.second << ")";
        segment_failure = stream.str();
        return std::vector<int>{};
      }

      int selected = candidates.front();
      std::vector<int> selected_connection;
      if (!faces.empty()) {
        const Point2d event_point{
          static_cast<double>(static_cast<long double>(start.first) +
                              events[interval] *
                                (static_cast<long double>(goal.first) - start.first)),
          static_cast<double>(static_cast<long double>(start.second) +
                              events[interval] *
                                (static_cast<long double>(goal.second) - start.second))};
        size_t best_size = std::numeric_limits<size_t>::max();
        for (const int candidate : candidates) {
          auto connection = connect_at_point(faces.back(), candidate, event_point);
          if (!connection.empty() && (connection.size() < best_size ||
                                      (connection.size() == best_size && candidate < selected))) {
            selected = candidate;
            best_size = connection.size();
            selected_connection = std::move(connection);
          }
        }
        if (selected_connection.empty()) {
          std::ostringstream stream;
          stream << "no incident free-face connection from triangle " << faces.back()
                 << " at parameter " << static_cast<double>(events[interval]);
          segment_failure = stream.str();
          return std::vector<int>{};
        }
        faces.insert(faces.end(), selected_connection.begin() + 1, selected_connection.end());
      }
      if (faces.empty() || faces.back() != selected)
        faces.push_back(selected);
    }
    return faces;
  };

  std::vector<int> raw_faces;
  for (size_t segment_index = 1; segment_index < path.size(); ++segment_index) {
    if (stop_token.poll())
      return stopped();
    auto current_faces = segment_faces(path[segment_index - 1], path[segment_index]);
    // segment_faces() and its connect_at_point() helper use an empty vector
    // for both geometric failure and cooperative stop.  Check the latched
    // token before interpreting that sentinel as a missing free-space sleeve.
    if (stop_token.poll())
      return stopped();
    if (current_faces.empty()) {
      return fail("Reference segment " + std::to_string(segment_index - 1) + "->" +
                  std::to_string(segment_index) + " leaves the free-triangle environment" +
                  (segment_failure.empty() ? std::string{} : ": " + segment_failure));
    }
    if (!raw_faces.empty()) {
      auto connection =
        connect_at_point(raw_faces.back(), current_faces.front(), path[segment_index - 1]);
      if (stop_token.poll())
        return stopped();
      if (connection.empty())
        return fail("Consecutive reference segments do not share a free triangle fan");
      raw_faces.insert(raw_faces.end(), connection.begin() + 1, connection.end());
    }
    if (raw_faces.empty()) {
      raw_faces.insert(raw_faces.end(), current_faces.begin(), current_faces.end());
    } else {
      const size_t offset = raw_faces.back() == current_faces.front() ? 1 : 0;
      raw_faces.insert(raw_faces.end(), current_faces.begin() + offset, current_faces.end());
    }
  }

  std::vector<int> reduced_faces;
  reduced_faces.reserve(raw_faces.size());
  for (const int face : raw_faces) {
    if (!reduced_faces.empty() && reduced_faces.back() == face)
      continue;
    if (reduced_faces.size() >= 2 && reduced_faces[reduced_faces.size() - 2] == face) {
      reduced_faces.pop_back();
      continue;
    }
    reduced_faces.push_back(face);
  }
  if (reduced_faces.empty())
    return fail("Reference path produced an empty triangle corridor");

  corridor.triangle_occurrences.reserve(reduced_faces.size());
  for (const int face : reduced_faces)
    corridor.triangle_occurrences.push_back(static_cast<uint32_t>(face));
  corridor.portals.reserve(reduced_faces.size() - 1);
  for (size_t index = 1; index < reduced_faces.size(); ++index) {
    const int from = reduced_faces[index - 1];
    const int to = reduced_faces[index];
    const auto& from_face = triangle_faces_[static_cast<size_t>(from)];
    int shared_edge = -1;
    for (size_t edge_index = 0; edge_index < 3; ++edge_index) {
      if (from_face.neighbors[edge_index] == to) {
        shared_edge = static_cast<int>(edge_index);
        break;
      }
    }
    if (shared_edge < 0)
      return fail("Reduced triangle corridor contains non-adjacent faces");
    const Point2d first = from_face.vertices[static_cast<size_t>(shared_edge)];
    const Point2d second = from_face.vertices[(static_cast<size_t>(shared_edge) + 1) % 3];
    DirectedPortal portal;
    portal.from_triangle = static_cast<uint32_t>(from);
    portal.to_triangle = static_cast<uint32_t>(to);
    // The common edge is stored in the boundary order of from_face. Unlike
    // a centroid-to-centroid test, this remains well defined for obtuse
    // adjacent triangles where both portal endpoints can lie on the same
    // side of the centroid line. For a counter-clockwise source face the
    // source interior is left of first->second, so crossing to the neighbor
    // makes second the portal's geometric left endpoint and first its right.
    // Reversing the traversal therefore swaps the endpoints exactly.
    const int face_orientation =
      exactOrientation(from_face.vertices[0], from_face.vertices[1], from_face.vertices[2]);
    if (face_orientation > 0) {
      portal.left = second;
      portal.right = first;
    } else {
      portal.left = first;
      portal.right = second;
    }
    corridor.portals.push_back(portal);
  }
  return OperationStatus::success;
}

HomotopyShorteningResult Polymap::shortenPathWithinHomotopy(const std::vector<Point2d>& reference,
                                                            const StopToken& stop_token) const {
  HomotopyShorteningResult result;
  const auto clear_uncertified_output = [&result]() {
    result.path.clear();
    result.path_cost = 0.0;
    result.collision_free = false;
    result.homotopy_preserved = false;
    result.locally_shortest = false;
  };
  if (stop_token.poll()) {
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled before tracing";
    return result;
  }
  std::string trace_error;
  const OperationStatus trace_status =
    traceFreeSpacePath(reference, result.corridor, stop_token, &trace_error);
  if (trace_status == OperationStatus::stopped) {
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled while tracing the reference";
    return result;
  }
  if (trace_status != OperationStatus::success) {
    result.status = reference.empty() ? HomotopyShorteningStatus::invalid_reference
                                      : HomotopyShorteningStatus::no_corridor;
    result.message = trace_error.empty() ? "Could not trace the reference path" : trace_error;
    return result;
  }
  std::string witness_error;
  if (!validateReducedDirectedPortalWitness(result.corridor, &witness_error)) {
    result.status = HomotopyShorteningStatus::failure;
    result.message = "Reference trace produced an invalid lifted-sleeve witness: " + witness_error;
    return result;
  }
  if (stop_token.poll()) {
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled before funnel execution";
    return result;
  }

  result.path = runFunnel(reference.front(), reference.back(), result.corridor.portals);
  if (result.path.empty() || !samePoint(result.path.front(), reference.front()) ||
      !samePoint(result.path.back(), reference.back())) {
    clear_uncertified_output();
    result.status = HomotopyShorteningStatus::failure;
    result.message = "Funnel result did not preserve the exact reference endpoints";
    return result;
  }
  result.path_cost = polylineLength(result.path);
  const double reference_cost = polylineLength(reference);
  const double cost_tolerance = 1.0e-10 * std::max({1.0, reference_cost, result.path_cost});
  if (result.path_cost > reference_cost + cost_tolerance) {
    std::ostringstream diagnostic;
    diagnostic << "Funnel result cost " << result.path_cost << " exceeds reference cost "
               << reference_cost << "; output=";
    for (const auto& point : result.path)
      diagnostic << "(" << point.first << "," << point.second << ")";
    diagnostic << "; portals=";
    for (const auto& portal : result.corridor.portals) {
      diagnostic << "[" << portal.from_triangle << "->" << portal.to_triangle << " L("
                 << portal.left.first << "," << portal.left.second << ") R(" << portal.right.first
                 << "," << portal.right.second << ")]";
    }
    clear_uncertified_output();
    result.status = HomotopyShorteningStatus::failure;
    result.message = diagnostic.str();
    return result;
  }

  std::string output_error;
  const OperationStatus output_status =
    traceFreeSpacePath(result.path, result.output_corridor, stop_token, &output_error);
  if (output_status == OperationStatus::stopped) {
    clear_uncertified_output();
    result.output_corridor = {};
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled during output validation";
    return result;
  }
  if (output_status != OperationStatus::success) {
    clear_uncertified_output();
    result.output_corridor = {};
    result.status = HomotopyShorteningStatus::failure;
    result.message = "Funnel result failed free-space validation: " + output_error;
    return result;
  }
  if (stop_token.poll()) {
    clear_uncertified_output();
    result.output_corridor = {};
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled before portal-witness certification";
    return result;
  }
  // The contour CDT adds no free-space Steiner vertices: every free face is
  // convex and every unconstrained shared edge is a contractible portal.
  // Contracting path fragments inside successive faces therefore maps a path
  // to its reduced ordered dual-edge walk.  Equality of the complete directed
  // portal-occurrence walk (not merely a set of face IDs) is a sufficient
  // relative-endpoint homotopy certificate, including repeated winding
  // cycles.  Re-tracing the emitted geometry makes this proof independent of
  // the funnel implementation itself.
  if (!sameReducedDirectedPortalWitness(result.corridor, result.output_corridor, &witness_error)) {
    clear_uncertified_output();
    result.status = HomotopyShorteningStatus::failure;
    result.message = "Funnel result changed the reduced lifted-sleeve witness: " + witness_error;
    return result;
  }
  if (stop_token.poll()) {
    clear_uncertified_output();
    result.output_corridor = {};
    result.status = HomotopyShorteningStatus::stopped;
    result.message = "Homotopy shortening was canceled after portal-witness certification";
    return result;
  }

  result.status = HomotopyShorteningStatus::success;
  result.collision_free = true;
  result.homotopy_preserved = true;
  result.locally_shortest = true;
  result.message =
    "Funnel completed and independently retraced the same reduced directed-portal witness";
  return result;
}

bool Polymap::isFacetInsideObstacle(int facet_idx) const {
  if (facet_idx < 0 || facet_idx >= static_cast<int>(facets_.size()))
    return false;
  const auto& f = facets_[facet_idx];
  double cx = (f[0].first + f[1].first + f[2].first) / 3.0;
  double cy = (f[0].second + f[1].second + f[2].second) / 3.0;
  for (size_t oi = 1; oi < obs_.size(); ++oi) {
    const auto& v = obs_[oi].ordered_vertices_;
    int n = static_cast<int>(v.size());
    bool inside = false;
    for (int k = 0, j = n - 1; k < n; j = k++) {
      if (((v[k].second > cy) != (v[j].second > cy)) &&
          (cx <
           (v[j].first - v[k].first) * (cy - v[k].second) / (v[j].second - v[k].second + 1e-15) +
             v[k].first))
        inside = !inside;
    }
    if (inside)
      return true;
  }
  return false;
}

std::vector<Polymap::CDTEdge> Polymap::getCDTEdges(size_t max_edges) const {
  std::vector<CDTEdge> edges;
  if (!cdt_ready_ || max_edges == 0)
    return edges;

  // Do not reserve the full triangulation: this method is also used by the
  // ROS visualization path, where the caller may have a deliberately small
  // marker budget.
  const size_t reserve_hint = std::min<size_t>(max_edges, 1024);
  edges.reserve(reserve_hint);
  for (auto eit = cdt_.finite_edges_begin(); eit != cdt_.finite_edges_end(); ++eit) {
    if (edges.size() >= max_edges)
      break;
    auto face = eit->first;
    int idx = eit->second;
    auto pa = face->vertex(cdt_.cw(idx))->point();
    auto pb = face->vertex(cdt_.ccw(idx))->point();
    CDTEdge e;
    e.a = {static_cast<int>(pa.x()), static_cast<int>(pa.y())};
    e.b = {static_cast<int>(pb.x()), static_cast<int>(pb.y())};
    e.is_constrained = cdt_.is_constrained(*eit);
    edges.push_back(e);
  }
  return edges;
}


}  // namespace raystar
