#pragma once

#include <rclcpp/rclcpp.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <nav_msgs/msg/path.hpp>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <visualization_msgs/msg/marker_array.hpp>
#include <raystar_interfaces/srv/get_raystar_paths.hpp>

#include <raystar/raystar_core.h>

namespace raystar
{

class RaystarNode : public rclcpp::Node
{
public:
  explicit RaystarNode(const rclcpp::NodeOptions& options = rclcpp::NodeOptions());

private:
  rclcpp::Service<raystar_interfaces::srv::GetRaystarPaths>::SharedPtr service_;

  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr non_homotopic_pub_;
  rclcpp::Publisher<visualization_msgs::msg::MarkerArray>::SharedPtr poly_obstacle_pub_;

  RaystarCore core_;

  void handleService(
    const std::shared_ptr<raystar_interfaces::srv::GetRaystarPaths::Request> request,
    std::shared_ptr<raystar_interfaces::srv::GetRaystarPaths::Response> response);

  GridMap occupancyGridToBinaryMap(
    const nav_msgs::msg::OccupancyGrid& grid, bool allow_unknown) const;

  void publishPolyObstacles(const Polymap& polymap, const GridMap& grid_map);
  void publishNonHomotopicPaths(
    const std::vector<PathSolution>& solutions, const GridMap& grid_map);

  nav_msgs::msg::Path buildPathMsg(
    const std::vector<std::pair<int, int>>& path, const GridMap& grid_map) const;
};

}  // namespace raystar
