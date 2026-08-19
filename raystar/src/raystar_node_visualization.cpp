#include "raystar_node.h"
#include "conservative_path_length.h"
#include "metric_bound_search.h"
#include "published_path_order.h"

#include <std_msgs/msg/color_rgba.hpp>
#include <raystar_interfaces/msg/debug_node.hpp>
#include <raystar_interfaces/msg/path_result.hpp>
#include <raystar_interfaces/msg/planning_result_info.hpp>
#include <rcl_interfaces/msg/integer_range.hpp>
#include <rcl_interfaces/msg/parameter_descriptor.hpp>
#include <rcl_interfaces/msg/set_parameters_result.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <type_traits>
#include <utility>
#include <exception>

#include "raystar_node_detail.h"

namespace raystar {

using namespace node_impl;

bool RaystarNode::buildPathMsg(const PathSolution& solution,
                               const GridMap& grid_map,
                               const std::string& frame_id,
                               size_t max_path_points,
                               nav_msgs::msg::Path& msg,
                               std::string& error) const {
  msg = nav_msgs::msg::Path{};
  error.clear();
  msg.header.stamp = now();
  msg.header.frame_id = frame_id;

  std::vector<Point2d> interpolated;
  if (!interpolateProjectedPath(solution, max_path_points, interpolated, error)) {
    return false;
  }
  msg.poses.reserve(interpolated.size());
  for (const auto& point : interpolated) {
    geometry_msgs::msg::PoseStamped pose;
    pose.header = msg.header;
    pose.pose.orientation.w = 1.0;
    double wx, wy;
    if (!continuousGridToWorld(grid_map, point, wx, wy)) {
      error = "path contains a point outside the finite world transform";
      msg.poses.clear();
      return false;
    }
    pose.pose.position.x = wx;
    pose.pose.position.y = wy;
    msg.poses.push_back(pose);
  }

  return true;
}

bool RaystarNode::buildTopologyPathMsg(const PathSolution& solution,
                                       const GridMap& grid_map,
                                       const std::string& frame_id,
                                       nav_msgs::msg::Path& msg,
                                       std::string& error) const {
  msg = nav_msgs::msg::Path{};
  error.clear();
  msg.header.stamp = now();
  msg.header.frame_id = frame_id;

  const auto projected = solution.projectedPath();
  if (projected.size() < 2) {
    error = "topology path must contain both endpoints";
    return false;
  }
  msg.poses.reserve(projected.size());
  for (const auto& point : projected) {
    geometry_msgs::msg::PoseStamped pose;
    pose.header = msg.header;
    pose.pose.orientation.w = 1.0;
    if (!continuousGridToWorld(grid_map, point, pose.pose.position.x, pose.pose.position.y)) {
      error = "topology path contains a point outside the finite world transform";
      msg.poses.clear();
      return false;
    }
    msg.poses.emplace_back(std::move(pose));
  }
  return true;
}

void RaystarNode::publishPolyObstacles(const Polymap& polymap,
                                       const GridMap& grid_map,
                                       const std::string& frame_id,
                                       size_t max_marker_bytes) {
  const auto stamp = now();
  auto array = makeMarkerSnapshot(frame_id, stamp);
  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = frame_id;
  marker.header.stamp = stamp;
  marker.ns = "polygons";
  marker.id = 0;
  marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
  marker.action = visualization_msgs::msg::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.scale.x = 0.05;
  marker.color.r = 1.0;
  marker.color.a = 1.0;

  const size_t topic_budget = markerTopicBudget(max_marker_bytes);
  const size_t point_limit = markerPointLimit(topic_budget);
  const size_t marker_limit = markerEntryLimit(topic_budget);
  size_t emitted_points = 0;
  size_t emitted_markers = 0;
  bool truncated = false;
  for (const auto& ob : polymap.obstacles()) {
    if (emitted_markers >= marker_limit) {
      truncated = true;
      break;
    }
    marker.id++;
    marker.points.clear();
    marker.colors.clear();

    for (auto it = ob.ordered_vertices_.begin(); it != ob.ordered_vertices_.end(); ++it) {
      auto nxt = std::next(it);
      if (nxt == ob.ordered_vertices_.end())
        nxt = ob.ordered_vertices_.begin();

      if (!canAppendCount(emitted_points, 2, point_limit)) {
        truncated = true;
        break;
      }

      double wx1, wy1, wx2, wy2;
      mapToWorld(grid_map,
                 static_cast<unsigned int>(it->first),
                 static_cast<unsigned int>(it->second),
                 wx1,
                 wy1);
      mapToWorld(grid_map,
                 static_cast<unsigned int>(nxt->first),
                 static_cast<unsigned int>(nxt->second),
                 wx2,
                 wy2);

      geometry_msgs::msg::Point p1, p2;
      p1.x = wx1;
      p1.y = wy1;
      p2.x = wx2;
      p2.y = wy2;
      marker.points.push_back(p1);
      marker.points.push_back(p2);

      std_msgs::msg::ColorRGBA c;
      c.r = 1.0;
      c.a = 1.0;
      marker.colors.push_back(c);
      marker.colors.push_back(c);
      emitted_points += 2;
    }
    if (!marker.points.empty()) {
      array.markers.push_back(marker);
      ++emitted_markers;
    }
    if (truncated)
      break;
  }
  (void)truncated;
  poly_obstacle_pub_->publish(array);
}

void RaystarNode::publishCDT(const Polymap& polymap,
                             const GridMap& grid_map,
                             const std::string& frame_id,
                             size_t max_marker_bytes) {
  const size_t topic_budget = markerTopicBudget(max_marker_bytes);
  const size_t point_limit = markerPointLimit(topic_budget);
  const size_t marker_limit = markerEntryLimit(topic_budget);
  const size_t edge_limit = point_limit / 2;
  auto cdt_edges = polymap.getCDTEdges(edge_limit);

  const auto stamp = now();
  auto array = makeMarkerSnapshot(frame_id, stamp);

  // Internal (non-constrained) edges — thin blue lines
  visualization_msgs::msg::Marker int_marker;
  int_marker.header.frame_id = frame_id;
  int_marker.header.stamp = stamp;
  int_marker.ns = "cdt_internal";
  int_marker.id = 0;
  int_marker.type = visualization_msgs::msg::Marker::LINE_LIST;
  int_marker.action = visualization_msgs::msg::Marker::ADD;
  int_marker.pose.orientation.w = 1.0;
  int_marker.scale.x = 0.015;
  int_marker.color.r = 0.3;
  int_marker.color.g = 0.5;
  int_marker.color.b = 0.8;
  int_marker.color.a = 0.4;

  // Constrained (obstacle) edges — thick red lines
  visualization_msgs::msg::Marker con_marker;
  con_marker.header.frame_id = frame_id;
  con_marker.header.stamp = stamp;
  con_marker.ns = "cdt_constrained";
  con_marker.id = 0;
  con_marker.type = visualization_msgs::msg::Marker::LINE_LIST;
  con_marker.action = visualization_msgs::msg::Marker::ADD;
  con_marker.pose.orientation.w = 1.0;
  con_marker.scale.x = 0.06;
  con_marker.color.r = 1.0;
  con_marker.color.g = 0.2;
  con_marker.color.b = 0.2;
  con_marker.color.a = 0.9;

  size_t emitted_points = 0;
  for (const auto& e : cdt_edges) {
    if (!canAppendCount(emitted_points, 2, point_limit))
      break;
    geometry_msgs::msg::Point p1, p2;
    // Promote integer vertices before multiplying by the ROS float
    // resolution.  Otherwise the multiplication is performed in float and
    // large grid indices can lose several low bits before being assigned to
    // the double-valued Marker point.
    const double resolution = static_cast<double>(grid_map.resolution);
    p1.x = grid_map.origin_x + static_cast<double>(e.a.first) * resolution;
    p1.y = grid_map.origin_y + static_cast<double>(e.a.second) * resolution;
    p2.x = grid_map.origin_x + static_cast<double>(e.b.first) * resolution;
    p2.y = grid_map.origin_y + static_cast<double>(e.b.second) * resolution;

    if (e.is_constrained) {
      con_marker.points.push_back(p1);
      con_marker.points.push_back(p2);
    } else {
      int_marker.points.push_back(p1);
      int_marker.points.push_back(p2);
    }
    emitted_points += 2;
  }

  if (marker_limit > 0)
    array.markers.push_back(int_marker);
  if (marker_limit > 1)
    array.markers.push_back(con_marker);
  cdt_pub_->publish(array);
}

void RaystarNode::publishNonHomotopicPaths(const std::vector<PathSolution>& solutions,
                                           const GridMap& grid_map,
                                           const std::string& frame_id,
                                           size_t max_path_points,
                                           size_t max_marker_bytes) {
  // This snapshot can be replayed long after planning.  A zero stamp asks
  // RViz for the latest transform instead of pinning the cached message to a
  // time that may have fallen out of a dynamic TF buffer.
  const auto stamp = rclcpp::Time(0, 0, get_clock()->get_clock_type());
  auto array = makeMarkerSnapshot(frame_id, stamp);
  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = frame_id;
  marker.header.stamp = stamp;
  marker.id = 0;
  marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
  marker.action = visualization_msgs::msg::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.scale.x = 0.1;
  marker.color.a = 1.0;

  const size_t topic_budget = markerTopicBudget(max_marker_bytes);
  const size_t point_limit = markerPointLimit(topic_budget);
  const size_t marker_limit = markerEntryLimit(topic_budget);
  size_t emitted_points = 0;
  size_t emitted_markers = 0;
  int num_div = static_cast<int>(std::ceil(std::sqrt(solutions.size())));
  int step = 100 / (num_div + 1);

  for (size_t si = 0; si < solutions.size(); ++si) {
    if (emitted_markers >= marker_limit || emitted_points >= point_limit)
      break;
    marker.id++;
    marker.ns = "path_" + std::to_string(si + 1);
    marker.points.clear();
    marker.colors.clear();

    int ri = (si / ((num_div + 1) * (num_div + 1)));
    int gi = (static_cast<int>(si) / (num_div + 1)) % (num_div + 1);
    int bi = static_cast<int>(si) % (num_div + 1);
    std_msgs::msg::ColorRGBA color;
    color.r = (100 + ri * step) / 255.0;
    color.g = (100 + gi * step) / 255.0;
    color.b = (100 + bi * step) / 255.0;
    color.a = 1.0;

    std::vector<Point2d> interpolated;
    std::string interpolation_error;
    const size_t remaining_points = point_limit - emitted_points;
    if (!interpolateProjectedPath(solutions[si],
                                  std::min(max_path_points, remaining_points),
                                  interpolated,
                                  interpolation_error)) {
      RCLCPP_WARN(get_logger(),
                  "Skipping path visualization because output limits rejected it: %s",
                  interpolation_error.c_str());
      continue;
    }
    for (const auto& point : interpolated) {
      double wx, wy;
      if (!continuousGridToWorld(grid_map, point, wx, wy))
        continue;
      geometry_msgs::msg::Point path_point;
      path_point.x = wx;
      path_point.y = wy;
      marker.points.push_back(path_point);
      marker.colors.push_back(color);
      ++emitted_points;
    }
    if (!marker.points.empty()) {
      array.markers.push_back(marker);
      ++emitted_markers;
    }
  }
  auto snapshot = std::make_shared<MarkerArray>(std::move(array));
  non_homotopic_pub_->publish(*snapshot);
  if (path_visualization_timer_)
    cached_path_visualization_ = std::move(snapshot);
}

void RaystarNode::clearVisualizationsLocked() noexcept {
  std::string frame_id;
  try {
    // Accepted request frames are bounded by kMaxFrameIdBytes.  If even this
    // small copy cannot be made, an empty frame is still safe for DELETEALL.
    frame_id = last_frame_id_;
  } catch (...) {
    frame_id.clear();
  }

  // Drop the previous search tree before the next request's map conversion.
  // clear() alone retains each Node/vector capacity and can leave a large
  // failed or invalid request resident indefinitely.
  core_.resetSearchState();
  cached_path_visualization_.reset();
  last_frame_id_.clear();

  // Publishing is normally non-throwing, but a middleware allocation failure
  // must not escape from an exception handler and terminate the executor.  A
  // failure on one topic must also not prevent best-effort clears on the
  // remaining topics.
  try {
    const auto clear_snapshot = makeMarkerSnapshot(frame_id, now());
    for (const auto& publisher :
         {non_homotopic_pub_, poly_obstacle_pub_, debug_tree_pub_, cdt_pub_}) {
      try {
        publisher->publish(clear_snapshot);
      } catch (...) {
        // Continue clearing the other durable topics.  The local caches and
        // Core state have already been released.
      }
    }
  } catch (...) {
    // Even the small DELETEALL snapshot could not be allocated.
  }
}

void RaystarNode::republishCachedPathVisualization() {
  std::unique_lock<std::mutex> planner_lock(planner_cache_mutex_, std::try_to_lock);
  if (!planner_lock.owns_lock() || planning_busy_.load(std::memory_order_acquire)) {
    return;
  }
  if (!cached_path_visualization_)
    return;

  // Keep publish inside planner_cache_mutex_.  Copying the shared_ptr and
  // publishing after unlock would allow a new request to publish DELETEALL
  // first and this old snapshot second, resurrecting stale paths.
  try {
    if (non_homotopic_pub_->get_subscription_count() == 0)
      return;
    non_homotopic_pub_->publish(*cached_path_visualization_);
  } catch (const std::exception& exception) {
    RCLCPP_WARN_THROTTLE(get_logger(),
                         *get_clock(),
                         5000,
                         "Could not republish cached path visualization: %s",
                         exception.what());
  } catch (...) {
    RCLCPP_WARN_THROTTLE(
      get_logger(), *get_clock(), 5000, "Could not republish cached path visualization");
  }
}

void RaystarNode::publishDebugTree(const std::vector<raystar::Node>& nodes,
                                   const GridMap& grid_map,
                                   const std::string& frame_id,
                                   size_t max_debug_nodes,
                                   size_t max_marker_bytes) {
  const auto stamp = now();
  auto array = makeMarkerSnapshot(frame_id, stamp);
  const size_t node_count = std::min(nodes.size(), max_debug_nodes);
  if (node_count == 0) {
    debug_tree_pub_->publish(array);
    return;
  }

  const size_t topic_budget = markerTopicBudget(max_marker_bytes);
  const size_t point_limit = markerPointLimit(topic_budget);
  const size_t marker_limit = markerEntryLimit(topic_budget);
  // Keep room for the aggregate edge and seed markers appended after the
  // per-node text/region markers.
  const size_t node_marker_limit = marker_limit > 2 ? marker_limit - 2 : 0;
  size_t emitted_points = 0;
  size_t emitted_markers = 0;

  double min_f = std::numeric_limits<double>::max();
  double max_f = std::numeric_limits<double>::lowest();
  const double resolution = static_cast<double>(grid_map.resolution);
  for (size_t i = 0; i < node_count; ++i) {
    const auto& n = nodes[i];
    const double f = (n.gCost() + n.hCost()) * resolution;
    if (f < min_f)
      min_f = f;
    if (f > max_f)
      max_f = f;
  }

  visualization_msgs::msg::Marker edge_marker;
  edge_marker.header.frame_id = frame_id;
  edge_marker.header.stamp = stamp;
  edge_marker.ns = "tree_edges";
  edge_marker.id = 0;
  edge_marker.type = visualization_msgs::msg::Marker::LINE_LIST;
  edge_marker.action = visualization_msgs::msg::Marker::ADD;
  edge_marker.pose.orientation.w = 1.0;
  edge_marker.scale.x = 0.02;
  edge_marker.color.r = 0.5;
  edge_marker.color.g = 0.5;
  edge_marker.color.b = 0.5;
  edge_marker.color.a = 0.6;

  visualization_msgs::msg::Marker seed_marker;
  seed_marker.header.frame_id = frame_id;
  seed_marker.header.stamp = stamp;
  seed_marker.ns = "seeds";
  seed_marker.id = 0;
  seed_marker.type = visualization_msgs::msg::Marker::POINTS;
  seed_marker.action = visualization_msgs::msg::Marker::ADD;
  seed_marker.pose.orientation.w = 1.0;
  seed_marker.scale.x = 0.12;
  seed_marker.scale.y = 0.12;

  visualization_msgs::msg::Marker text_marker;
  text_marker.header.frame_id = frame_id;
  text_marker.header.stamp = stamp;
  text_marker.ns = "costs";
  text_marker.type = visualization_msgs::msg::Marker::TEXT_VIEW_FACING;
  text_marker.action = visualization_msgs::msg::Marker::ADD;
  text_marker.pose.orientation.w = 1.0;
  text_marker.scale.z = 0.1;
  text_marker.color.r = 1.0;
  text_marker.color.g = 1.0;
  text_marker.color.b = 1.0;
  text_marker.color.a = 1.0;

  for (size_t i = 0; i < node_count; ++i) {
    if (emitted_markers >= node_marker_limit || emitted_points >= point_limit)
      break;
    const auto& n = nodes[i];
    double wx, wy;
    mapToWorld(grid_map,
               static_cast<unsigned int>(n.seed().first),
               static_cast<unsigned int>(n.seed().second),
               wx,
               wy);

    const double f = (n.gCost() + n.hCost()) * resolution;
    double t = (max_f > min_f) ? (f - min_f) / (max_f - min_f) : 0.0;

    std_msgs::msg::ColorRGBA color;
    color.r = t;
    color.g = 1.0 - t;
    color.b = 0.2;
    color.a = 1.0;

    geometry_msgs::msg::Point p;
    p.x = wx;
    p.y = wy;
    if (!canAppendCount(emitted_points, 1, point_limit))
      break;
    seed_marker.points.push_back(p);
    seed_marker.colors.push_back(color);
    ++emitted_points;

    if (n.parentIndex() >= 0 && static_cast<size_t>(n.parentIndex()) < node_count) {
      const auto& parent = nodes[static_cast<size_t>(n.parentIndex())];
      double pwx, pwy;
      mapToWorld(grid_map,
                 static_cast<unsigned int>(parent.seed().first),
                 static_cast<unsigned int>(parent.seed().second),
                 pwx,
                 pwy);
      geometry_msgs::msg::Point pp;
      pp.x = pwx;
      pp.y = pwy;
      if (canAppendCount(emitted_points, 2, point_limit)) {
        edge_marker.points.push_back(pp);
        edge_marker.points.push_back(p);
        emitted_points += 2;
      }
    }

    text_marker.id = static_cast<int>(i);
    text_marker.pose.position.x = wx;
    text_marker.pose.position.y = wy + 0.2;
    text_marker.pose.position.z = 0.0;
    std::ostringstream oss;
    oss << "N" << n.index() << " G=" << std::fixed << std::setprecision(1) << n.gCost() * resolution
        << " F=" << std::fixed << std::setprecision(1) << f;
    text_marker.text = oss.str();
    if (emitted_markers >= node_marker_limit)
      break;
    array.markers.push_back(text_marker);
    ++emitted_markers;

    if (!n.visibility().empty() && emitted_markers < node_marker_limit &&
        emitted_points < point_limit) {
      visualization_msgs::msg::Marker visreg_marker;
      visreg_marker.header.frame_id = frame_id;
      visreg_marker.header.stamp = stamp;
      visreg_marker.ns = "node_" + std::to_string(n.index());
      visreg_marker.id = 0;
      visreg_marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
      visreg_marker.action = visualization_msgs::msg::Marker::ADD;
      visreg_marker.pose.orientation.w = 1.0;
      visreg_marker.scale.x = 0.03;
      visreg_marker.color.g = 0.8;
      visreg_marker.color.a = 0.5;
      for (const auto& v : n.visibility()) {
        if (!canAppendCount(emitted_points, 1, point_limit))
          break;
        geometry_msgs::msg::Point vp;
        vp.x = grid_map.origin_x + v.first * grid_map.resolution;
        vp.y = grid_map.origin_y + v.second * grid_map.resolution;
        visreg_marker.points.push_back(vp);
        ++emitted_points;
      }
      if (n.index() > 0) {
        if (canAppendCount(emitted_points, 1, point_limit)) {
          geometry_msgs::msg::Point sp;
          sp.x = grid_map.origin_x + static_cast<double>(n.seed().first) * grid_map.resolution;
          sp.y = grid_map.origin_y + static_cast<double>(n.seed().second) * grid_map.resolution;
          visreg_marker.points.push_back(sp);
          ++emitted_points;
        }
      }
      if (visreg_marker.points.size() > 2 && canAppendCount(emitted_points, 1, point_limit)) {
        visreg_marker.points.push_back(visreg_marker.points.front());
        ++emitted_points;
      }
      if (!visreg_marker.points.empty() && emitted_markers < node_marker_limit) {
        array.markers.push_back(visreg_marker);
        ++emitted_markers;
      }

      if (!n.fullVisibility().empty() && n.fullVisibility().size() != n.visibility().size() &&
          emitted_markers < node_marker_limit && emitted_points < point_limit) {
        visualization_msgs::msg::Marker full_visreg_marker;
        full_visreg_marker.header.frame_id = frame_id;
        full_visreg_marker.header.stamp = stamp;
        full_visreg_marker.ns = "full_visreg_" + std::to_string(n.index());
        full_visreg_marker.id = 0;
        full_visreg_marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
        full_visreg_marker.action = visualization_msgs::msg::Marker::ADD;
        full_visreg_marker.pose.orientation.w = 1.0;
        full_visreg_marker.scale.x = 0.02;
        full_visreg_marker.color.r = 0.8;
        full_visreg_marker.color.g = 0.6;
        full_visreg_marker.color.b = 0.2;
        full_visreg_marker.color.a = 0.35;
        for (const auto& v : n.fullVisibility()) {
          if (!canAppendCount(emitted_points, 1, point_limit))
            break;
          geometry_msgs::msg::Point vp;
          vp.x = grid_map.origin_x + v.first * grid_map.resolution;
          vp.y = grid_map.origin_y + v.second * grid_map.resolution;
          full_visreg_marker.points.push_back(vp);
          ++emitted_points;
        }
        if (n.index() > 0) {
          if (canAppendCount(emitted_points, 1, point_limit)) {
            geometry_msgs::msg::Point sp;
            sp.x = grid_map.origin_x + static_cast<double>(n.seed().first) * grid_map.resolution;
            sp.y = grid_map.origin_y + static_cast<double>(n.seed().second) * grid_map.resolution;
            full_visreg_marker.points.push_back(sp);
            ++emitted_points;
          }
        }
        if (full_visreg_marker.points.size() > 2 &&
            canAppendCount(emitted_points, 1, point_limit)) {
          full_visreg_marker.points.push_back(full_visreg_marker.points.front());
          ++emitted_points;
        }
        if (!full_visreg_marker.points.empty() && emitted_markers < node_marker_limit) {
          array.markers.push_back(full_visreg_marker);
          ++emitted_markers;
        }
      }

      if (n.localShortestPath().size() >= 2 && emitted_markers < node_marker_limit &&
          emitted_points < point_limit) {
        visualization_msgs::msg::Marker tpath_marker;
        tpath_marker.header.frame_id = frame_id;
        tpath_marker.header.stamp = stamp;
        tpath_marker.ns = "node_" + std::to_string(n.index());
        tpath_marker.id = 1;
        tpath_marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
        tpath_marker.action = visualization_msgs::msg::Marker::ADD;
        tpath_marker.pose.orientation.w = 1.0;
        tpath_marker.scale.x = 0.04;
        tpath_marker.color.r = 0.0;
        tpath_marker.color.g = 0.5;
        tpath_marker.color.b = 1.0;
        tpath_marker.color.a = 0.8;
        for (const auto& wp : n.localShortestPath()) {
          if (!canAppendCount(emitted_points, 1, point_limit))
            break;
          double tpx, tpy;
          mapToWorld(grid_map,
                     static_cast<unsigned int>(wp.first),
                     static_cast<unsigned int>(wp.second),
                     tpx,
                     tpy);
          geometry_msgs::msg::Point tp;
          tp.x = tpx;
          tp.y = tpy;
          tpath_marker.points.push_back(tp);
          ++emitted_points;
        }
        if (!tpath_marker.points.empty() && emitted_markers < node_marker_limit) {
          array.markers.push_back(tpath_marker);
          ++emitted_markers;
        }
      }
    }
  }

  if (emitted_markers < marker_limit) {
    array.markers.push_back(edge_marker);
    ++emitted_markers;
  }
  if (emitted_markers < marker_limit)
    array.markers.push_back(seed_marker);

  debug_tree_pub_->publish(array);
}


}  // namespace raystar
