#include "raystar_node.h"

#include <cmath>
#include <chrono>

namespace raystar
{

RaystarNode::RaystarNode(const rclcpp::NodeOptions& options)
  : Node("raystar", options)
{
  service_ = create_service<raystar_interfaces::srv::GetRaystarPaths>(
    "~/get_raystar_paths",
    std::bind(&RaystarNode::handleService, this,
      std::placeholders::_1, std::placeholders::_2));

  non_homotopic_pub_ = create_publisher<visualization_msgs::msg::MarkerArray>(
    "~/non_homotopic_paths", 1);
  poly_obstacle_pub_ = create_publisher<visualization_msgs::msg::MarkerArray>(
    "~/poly_obstacles", 1);

  RCLCPP_INFO(get_logger(), "Raystar node initialized");
}

GridMap RaystarNode::occupancyGridToBinaryMap(
  const nav_msgs::msg::OccupancyGrid& grid, bool allow_unknown) const
{
  GridMap map;
  map.width = grid.info.width;
  map.height = grid.info.height;
  map.resolution = grid.info.resolution;
  map.origin_x = grid.info.origin.position.x;
  map.origin_y = grid.info.origin.position.y;

  map.data.resize(map.width * map.height);
  for (size_t i = 0; i < grid.data.size(); ++i)
  {
    int8_t val = grid.data[i];
    if (val < 0) {
      map.data[i] = allow_unknown ? 0 : 1;
    } else if (val > 0) {
      map.data[i] = 1;
    } else {
      map.data[i] = 0;
    }
  }
  return map;
}

void RaystarNode::handleService(
  const std::shared_ptr<raystar_interfaces::srv::GetRaystarPaths::Request> request,
  std::shared_ptr<raystar_interfaces::srv::GetRaystarPaths::Response> response)
{
  auto& grid = request->map;

  if (grid.data.empty() || grid.info.width == 0 || grid.info.height == 0) {
    response->success = false;
    response->message = "Invalid map: empty or zero dimensions";
    return;
  }

  GridMap work_map = occupancyGridToBinaryMap(grid, request->allow_unknown);

  unsigned int start_mx, start_my, goal_mx, goal_my;
  if (!worldToMap(work_map,
      request->start.pose.position.x, request->start.pose.position.y,
      start_mx, start_my)) {
    response->success = false;
    response->message = "Start position is outside the map";
    return;
  }
  if (!worldToMap(work_map,
      request->goal.pose.position.x, request->goal.pose.position.y,
      goal_mx, goal_my)) {
    response->success = false;
    response->message = "Goal position is outside the map";
    return;
  }

  int start_x = static_cast<int>(start_mx);
  int start_y = static_cast<int>(start_my);
  int goal_x = static_cast<int>(goal_mx);
  int goal_y = static_cast<int>(goal_my);

  RCLCPP_INFO(get_logger(),
    "Planning: start=(%d,%d) goal=(%d,%d) K=%d allow_self_crossing=%s",
    start_x, start_y, goal_x, goal_y, request->k,
    request->allow_self_crossing ? "true" : "false");

  auto result = core_.plan(work_map, start_x, start_y, goal_x, goal_y,
    request->k, request->allow_self_crossing);

  RCLCPP_INFO(get_logger(),
    "Planning done: map=%.1fms plan=%.1fms paths=%zu",
    result.map_time_ms, result.plan_time_ms, result.path_solutions.size());

  response->success = result.success;
  response->message = result.message;

  for (auto& sol : result.path_solutions)
  {
    response->paths.push_back(buildPathMsg(sol.path_, work_map));
    response->path_costs.push_back(sol.path_cost_);
  }

  if (result.polymap) {
    publishPolyObstacles(*result.polymap, work_map);
    publishNonHomotopicPaths(result.path_solutions, work_map);
  }
}

nav_msgs::msg::Path RaystarNode::buildPathMsg(
  const std::vector<std::pair<int, int>>& path, const GridMap& grid_map) const
{
  nav_msgs::msg::Path msg;
  msg.header.stamp = now();
  msg.header.frame_id = "map";

  std::vector<std::pair<double, double>> interpolated;
  if (path.size() >= 2)
  {
    for (size_t i = 0; i < path.size() - 1; ++i)
    {
      auto& prev = path[i];
      auto& next = path[i + 1];
      double dis = std::hypot(prev.first - next.first, prev.second - next.second);
      int count = static_cast<int>(std::ceil(dis));
      for (int j = 0; j < count; ++j)
      {
        interpolated.emplace_back(
          (prev.first * (count - j) + next.first * j) / static_cast<double>(count),
          (prev.second * (count - j) + next.second * j) / static_cast<double>(count));
      }
    }
    interpolated.push_back({static_cast<double>(path.back().first),
      static_cast<double>(path.back().second)});
  }
  else if (path.size() == 1)
  {
    interpolated.push_back({static_cast<double>(path[0].first),
      static_cast<double>(path[0].second)});
  }

  for (auto& pt : interpolated)
  {
    geometry_msgs::msg::PoseStamped pose;
    pose.header = msg.header;
    pose.pose.orientation.w = 1.0;
    double wx, wy;
    mapToWorld(grid_map,
      static_cast<unsigned int>(std::round(pt.first)),
      static_cast<unsigned int>(std::round(pt.second)),
      wx, wy);
    pose.pose.position.x = wx;
    pose.pose.position.y = wy;
    msg.poses.push_back(pose);
  }

  return msg;
}

void RaystarNode::publishPolyObstacles(const Polymap& polymap, const GridMap& grid_map)
{
  if (poly_obstacle_pub_->get_subscription_count() == 0) return;

  visualization_msgs::msg::MarkerArray array;
  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = "map";
  marker.header.stamp = now();
  marker.ns = "polygons";
  marker.id = 0;
  marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
  marker.action = visualization_msgs::msg::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.scale.x = 0.1;
  marker.color.r = 1.0;
  marker.color.a = 1.0;

  for (auto& ob : polymap.obs_)
  {
    marker.id++;
    marker.points.clear();
    marker.colors.clear();

    for (auto it = ob.detailed_ordered_vertices_.begin();
      it != ob.detailed_ordered_vertices_.end(); ++it)
    {
      auto nxt = std::next(it);
      if (nxt == ob.detailed_ordered_vertices_.end())
        nxt = ob.detailed_ordered_vertices_.begin();

      double wx1, wy1, wx2, wy2;
      mapToWorld(grid_map,
        static_cast<unsigned int>(it->first),
        static_cast<unsigned int>(it->second), wx1, wy1);
      mapToWorld(grid_map,
        static_cast<unsigned int>(nxt->first),
        static_cast<unsigned int>(nxt->second), wx2, wy2);

      geometry_msgs::msg::Point p1, p2;
      p1.x = wx1; p1.y = wy1;
      p2.x = wx2; p2.y = wy2;
      marker.points.push_back(p1);
      marker.points.push_back(p2);

      std_msgs::msg::ColorRGBA c;
      c.r = 1.0; c.a = 1.0;
      marker.colors.push_back(c);
      marker.colors.push_back(c);
    }
    array.markers.push_back(marker);
  }
  poly_obstacle_pub_->publish(array);
}

void RaystarNode::publishNonHomotopicPaths(
  const std::vector<PathSolution>& solutions, const GridMap& grid_map)
{
  if (non_homotopic_pub_->get_subscription_count() == 0) return;

  visualization_msgs::msg::MarkerArray array;
  visualization_msgs::msg::Marker marker;
  marker.header.frame_id = "map";
  marker.header.stamp = now();
  marker.ns = "non_homotopic_paths";
  marker.id = 0;
  marker.type = visualization_msgs::msg::Marker::LINE_STRIP;
  marker.action = visualization_msgs::msg::Marker::ADD;
  marker.pose.orientation.w = 1.0;
  marker.scale.x = 0.1;
  marker.color.a = 1.0;

  int num_div = static_cast<int>(std::ceil(std::sqrt(solutions.size())));
  int step = 100 / (num_div + 1);

  for (size_t si = 0; si < solutions.size(); ++si)
  {
    marker.id++;
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

    auto& path = solutions[si].path_;
    for (auto it = path.begin(); it != path.end() - 1; ++it)
    {
      auto nxt = std::next(it);
      double wx1, wy1, wx2, wy2;
      mapToWorld(grid_map,
        static_cast<unsigned int>(it->first),
        static_cast<unsigned int>(it->second), wx1, wy1);
      mapToWorld(grid_map,
        static_cast<unsigned int>(nxt->first),
        static_cast<unsigned int>(nxt->second), wx2, wy2);

      geometry_msgs::msg::Point p1, p2;
      p1.x = wx1; p1.y = wy1;
      p2.x = wx2; p2.y = wy2;
      marker.points.push_back(p1);
      marker.points.push_back(p2);
      marker.colors.push_back(color);
      marker.colors.push_back(color);
    }
    array.markers.push_back(marker);
  }
  non_homotopic_pub_->publish(array);
}

}  // namespace raystar
