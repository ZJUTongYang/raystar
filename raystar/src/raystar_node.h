#pragma once

#include <rclcpp/rclcpp.hpp>
#include <rclcpp_action/rclcpp_action.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <nav_msgs/msg/path.hpp>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <visualization_msgs/msg/marker_array.hpp>
#include <raystar_interfaces/action/plan_raystar_paths.hpp>
#include <raystar_interfaces/map_identity.hpp>
#include <raystar_interfaces/msg/map_status.hpp>
#include <raystar_interfaces/srv/get_raystar_paths.hpp>

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>

#include <raystar/raystar_core.h>

namespace raystar
{

class RaystarNode : public rclcpp::Node
{
public:
  explicit RaystarNode(const rclcpp::NodeOptions& options = rclcpp::NodeOptions());
  ~RaystarNode() override;

private:
  using PlanAction = raystar_interfaces::action::PlanRaystarPaths;
  using PlanGoalHandle = rclcpp_action::ServerGoalHandle<PlanAction>;
  using MarkerArray = visualization_msgs::msg::MarkerArray;

  // Keep one immutable per-request configuration snapshot.  The same
  // PlanningLimits object is passed to Core and also drives ROS response/debug
  // admission, while the map threshold is used by the one binary conversion
  // for that request.
  struct RequestConfiguration
  {
    PlanningLimits planning;
    int occupied_threshold = 0;
  };

  // The callback is intentionally validation-only.  Keeping the handle alive
  // registers it for the node lifetime; it must not acquire planner locks or
  // mutate configuration while rclcpp holds its parameter mutex.
  rclcpp::node_interfaces::OnSetParametersCallbackHandle::SharedPtr
    parameter_validation_callback_;

  rclcpp::Service<raystar_interfaces::srv::GetRaystarPaths>::SharedPtr service_;
  rclcpp_action::Server<PlanAction>::SharedPtr action_server_;
  rclcpp::Subscription<nav_msgs::msg::OccupancyGrid>::SharedPtr map_subscription_;
  rclcpp::Publisher<raystar_interfaces::msg::MapStatus>::SharedPtr
    map_status_publisher_;

  rclcpp::Publisher<MarkerArray>::SharedPtr non_homotopic_pub_;
  rclcpp::Publisher<MarkerArray>::SharedPtr poly_obstacle_pub_;
  rclcpp::Publisher<MarkerArray>::SharedPtr debug_tree_pub_;
  rclcpp::Publisher<MarkerArray>::SharedPtr cdt_pub_;

  rclcpp::TimerBase::SharedPtr path_visualization_timer_;

  RaystarCore core_;

  // Service and Action share the same stateful Core instance.  Admission is
  // deliberately capacity one: a second request is rejected instead of
  // waiting in an unbounded executor/worker queue.
  std::atomic<bool> planning_busy_{false};
  std::atomic<bool> shutting_down_{false};
  std::mutex planner_cache_mutex_;

  struct CachedMap
  {
    nav_msgs::msg::OccupancyGrid::ConstSharedPtr message;
    raystar_interfaces::MapId id;
    std::uint64_t generation{0};
  };
  mutable std::mutex map_cache_mutex_;
  CachedMap cached_map_;

  // Action callbacks remain non-blocking.  The one admitted Action goal runs
  // on this node-owned worker, leaving the executor free to process cancel.
  struct ActionJob
  {
    std::shared_ptr<PlanGoalHandle> goal_handle;
    std::shared_ptr<std::atomic<bool>> cancel_requested;
  };

  std::thread action_worker_;
  std::mutex action_worker_mutex_;
  std::condition_variable action_worker_cv_;
  std::optional<ActionJob> pending_action_job_;
  bool stop_action_worker_{false};
  std::mutex action_state_mutex_;
  rclcpp_action::GoalUUID active_goal_id_{};
  bool active_goal_reserved_{false};
  std::shared_ptr<std::atomic<bool>> active_goal_cancel_;

  // RViz Humble deletes a disabled Marker namespace and does not request a
  // transient-local replay when that namespace is enabled again.  Retain the
  // already-built, bounded path snapshot so the timer can republish it without
  // repeating interpolation or Marker construction.
  std::shared_ptr<const MarkerArray> cached_path_visualization_;
  std::string last_frame_id_;

  void handleService(
    const std::shared_ptr<raystar_interfaces::srv::GetRaystarPaths::Request> request,
    std::shared_ptr<raystar_interfaces::srv::GetRaystarPaths::Response> response);

  rclcpp_action::GoalResponse handleActionGoal(
    const rclcpp_action::GoalUUID& uuid,
    std::shared_ptr<const PlanAction::Goal> goal);
  rclcpp_action::CancelResponse handleActionCancel(
    const std::shared_ptr<PlanGoalHandle> goal_handle);
  void handleActionAccepted(
    const std::shared_ptr<PlanGoalHandle> goal_handle);
  void executeAction(
    const std::shared_ptr<PlanGoalHandle> goal_handle,
    const std::shared_ptr<std::atomic<bool>>& cancel_requested) noexcept;
  void actionWorkerLoop() noexcept;

  template<typename RequestT, typename ResponseT>
  void executePlanning(
    const RequestT& request, ResponseT& response,
    const nav_msgs::msg::OccupancyGrid& grid,
    const raystar_interfaces::MapId& map_id,
    const StopPredicate& stop_requested) noexcept;

  void handleMap(
    nav_msgs::msg::OccupancyGrid::ConstSharedPtr map);
  bool resolveCachedMap(
    const raystar_interfaces::MapId& requested_id,
    nav_msgs::msg::OccupancyGrid::ConstSharedPtr& map,
    std::string& error) const;

  bool occupancyGridToBinaryMap(
    const nav_msgs::msg::OccupancyGrid& grid, bool allow_unknown,
    const RequestConfiguration& configuration,
    const StopPredicate& stop_requested,
    GridMap& output, std::string& error) const;
  bool loadRequestConfiguration(
    RequestConfiguration& configuration, std::string& error) const;

  void publishPolyObstacles(
    const Polymap& polymap, const GridMap& grid_map,
    const std::string& frame_id, size_t max_marker_bytes);
  void publishCDT(
    const Polymap& polymap, const GridMap& grid_map,
    const std::string& frame_id, size_t max_marker_bytes);
  void publishNonHomotopicPaths(
    const std::vector<PathSolution>& solutions, const GridMap& grid_map,
    const std::string& frame_id, size_t max_path_points,
    size_t max_marker_bytes);
  void publishDebugTree(
    const std::vector<raystar::Node>& nodes, const GridMap& grid_map,
    const std::string& frame_id, size_t max_debug_nodes,
    size_t max_marker_bytes);

  bool buildPathMsg(
    const PathSolution& solution, const GridMap& grid_map,
    const std::string& frame_id, size_t max_path_points,
    nav_msgs::msg::Path& message, std::string& error) const;

  // planner_cache_mutex_ must be held by the caller.
  void clearVisualizationsLocked() noexcept;
  void republishCachedPathVisualization();
};

}  // namespace raystar
