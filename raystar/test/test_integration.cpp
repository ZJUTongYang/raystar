#include <gtest/gtest.h>
#include <rclcpp/rclcpp.hpp>
#include <rclcpp_action/rclcpp_action.hpp>
#include <rclcpp/parameter_client.hpp>
#include <rcl_interfaces/msg/parameter_descriptor.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <visualization_msgs/msg/marker_array.hpp>
#include <raystar_interfaces/action/plan_raystar_paths.hpp>
#include <raystar_interfaces/map_identity.hpp>
#include <raystar_interfaces/msg/map_status.hpp>
#include <raystar_interfaces/msg/planning_result_info.hpp>
#include <raystar_interfaces/srv/get_raystar_paths.hpp>

#include "path_collision_oracle.h"
#include "small_grid_property.h"

#include <sys/wait.h>
#include <unistd.h>
#include <fcntl.h>

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <functional>
#include <limits>
#include <string>
#include <system_error>
#include <thread>
#include <vector>

using namespace std::chrono_literals;
using RaystarService = raystar_interfaces::srv::GetRaystarPaths;
using RaystarAction = raystar_interfaces::action::PlanRaystarPaths;
using RaystarGoalHandle = rclcpp_action::ClientGoalHandle<RaystarAction>;
using PlanningResultInfo = raystar_interfaces::msg::PlanningResultInfo;

static nav_msgs::msg::OccupancyGrid makeTestGrid()
{
  nav_msgs::msg::OccupancyGrid grid;
  grid.header.frame_id = "map";
  grid.info.width = 30;
  grid.info.height = 30;
  grid.info.resolution = 1.0f;
  grid.info.origin.position.x = 0.0;
  grid.info.origin.position.y = 0.0;
  grid.info.origin.orientation.w = 1.0;
  grid.data.resize(900, 0);

  for (unsigned int y = 10; y < 20; ++y)
    for (unsigned int x = 10; x < 20; ++x)
      grid.data[y * 30 + x] = 100;

  return grid;
}

static RaystarService::Request::SharedPtr makeTestRequest()
{
  auto request = std::make_shared<RaystarService::Request>();
  request->map = makeTestGrid();
  request->start.header.frame_id = "map";
  request->start.pose.position.x = 5.5;
  request->start.pose.position.y = 15.5;
  request->start.pose.orientation.w = 1.0;
  request->goal.header.frame_id = "map";
  request->goal.pose.position.x = 25.5;
  request->goal.pose.position.y = 15.5;
  request->goal.pose.orientation.w = 1.0;
  request->k = 3;
  request->allow_self_crossing = false;
  request->allow_unknown = false;
  request->include_debug = false;
  return request;
}

static RaystarService::Request::SharedPtr makeVerticalBarrierRequest(
  int8_t occupancy_value, bool allow_unknown = false)
{
  auto request = makeTestRequest();
  std::fill(request->map.data.begin(), request->map.data.end(), 0);
  const size_t width = static_cast<size_t>(request->map.info.width);
  const size_t height = static_cast<size_t>(request->map.info.height);
  const size_t barrier_x = width / 2;
  for (size_t y = 0; y < height; ++y)
    request->map.data[y * width + barrier_x] = occupancy_value;
  request->k = 1;
  request->allow_unknown = allow_unknown;
  return request;
}

static RaystarAction::Goal makeTestActionGoal()
{
  auto request = makeTestRequest();
  RaystarAction::Goal goal;
  goal.map_id = raystar_interfaces::computeMapId(request->map);
  goal.start = std::move(request->start);
  goal.goal = std::move(request->goal);
  goal.k = request->k;
  goal.allow_self_crossing = request->allow_self_crossing;
  goal.allow_unknown = request->allow_unknown;
  goal.include_debug = false;
  return goal;
}

static nav_msgs::msg::OccupancyGrid makeLongRunningActionMap()
{
  constexpr unsigned int width = 768;
  constexpr unsigned int height = 768;

  auto map = makeTestGrid();
  map.info.width = width;
  map.info.height = height;
  map.data.assign(static_cast<size_t>(width) * height, 0);

  // Thousands of isolated one-cell obstacles make contour/CDT construction
  // decisively longer than the Action goal/cancel round trips, while keeping
  // one connected free-space component and avoiding diagonal vertex contact.
  for (unsigned int y = 4; y + 4 < height; y += 4)
  {
    for (unsigned int x = 4; x + 4 < width; x += 4)
      map.data[static_cast<size_t>(y) * width + x] = 100;
  }
  return map;
}

static RaystarAction::Goal makeLongRunningActionGoal()
{
  constexpr unsigned int width = 768;
  constexpr unsigned int height = 768;
  const auto map = makeLongRunningActionMap();
  auto goal = makeTestActionGoal();
  goal.map_id = raystar_interfaces::computeMapId(map);
  goal.start.pose.position.x = 2.5;
  goal.start.pose.position.y = 2.5;
  goal.goal.pose.position.x = static_cast<double>(width) - 2.5;
  goal.goal.pose.position.y = static_cast<double>(height) - 2.5;
  goal.k = 100;
  return goal;
}

template<typename FutureT>
static bool waitForFuture(
  rclcpp::executors::SingleThreadedExecutor& executor,
  FutureT& future, std::chrono::milliseconds timeout)
{
  const auto deadline = std::chrono::steady_clock::now() + timeout;
  while (std::chrono::steady_clock::now() < deadline)
  {
    if (future.wait_for(0s) == std::future_status::ready)
      return true;
    executor.spin_some(20ms);
  }
  executor.spin_some(20ms);
  return future.wait_for(0s) == std::future_status::ready;
}

static bool cacheMapAndWait(
  rclcpp::executors::SingleThreadedExecutor& executor,
  const rclcpp::Node::SharedPtr& node,
  const nav_msgs::msg::OccupancyGrid& map,
  const std::string& server_namespace,
  std::chrono::milliseconds timeout = 5s)
{
  const auto expected_id = raystar_interfaces::computeMapId(map);
  const auto status_topic = server_namespace + "/raystar/map_status";
  const auto map_topic = server_namespace + "/map";
  bool admitted = false;
  auto status_subscription =
    node->create_subscription<raystar_interfaces::msg::MapStatus>(
      status_topic, rclcpp::QoS(rclcpp::KeepLast(1)).transient_local(),
      [&](raystar_interfaces::msg::MapStatus::ConstSharedPtr status) {
        admitted = status &&
          raystar_interfaces::mapIdsEqual(status->map_id, expected_id);
      });
  auto map_publisher =
    node->create_publisher<nav_msgs::msg::OccupancyGrid>(
      map_topic, rclcpp::QoS(rclcpp::KeepLast(1)).transient_local());

  const auto deadline = std::chrono::steady_clock::now() + timeout;
  auto next_publish = std::chrono::steady_clock::time_point::min();
  while (!admitted && std::chrono::steady_clock::now() < deadline)
  {
    const auto now = std::chrono::steady_clock::now();
    if (now >= next_publish)
    {
      map_publisher->publish(map);
      next_publish = now + 100ms;
    }
    executor.spin_some(20ms);
    std::this_thread::sleep_for(1ms);
  }
  (void)status_subscription;
  return admitted;
}

// Every integration-test invocation gets a private ROS namespace.  This is
// intentionally based on the parent PID rather than a fixed literal so that
// two ctest jobs can run concurrently without sharing services, actions or
// cached-map topics.  Child nodes inherit the same namespace prefix.
static const std::string& integrationNamespacePrefix()
{
  static const std::string prefix =
    "/raystar_integration_" + std::to_string(static_cast<long long>(getpid()));
  return prefix;
}

static std::string integrationNamespace(const std::string& scenario)
{
  return integrationNamespacePrefix() + "/" + scenario;
}

static std::string raystarEndpoint(
  const std::string& server_namespace, const std::string& suffix = {})
{
  return server_namespace + "/raystar" +
    (suffix.empty() ? std::string{} : "/" + suffix);
}

static std::string mainServerNamespace()
{
  return integrationNamespace("main");
}

struct SpawnedNode
{
  pid_t pid{-1};
  std::string node_namespace;
  int launch_errno{0};
};

// Start raystar_node and synchronously distinguish a successful exec from an
// immediate child-side failure.  The close-on-exec pipe reaches EOF only after
// exec has replaced the child image; otherwise the child writes errno before
// exiting.  This avoids reporting a discovery timeout when the executable (or
// its dynamic loader) was actually missing.
static SpawnedNode spawnRaystarNode(
  const std::string& scenario,
  const std::vector<std::string>& parameter_overrides = {})
{
  SpawnedNode spawned;
  spawned.node_namespace = integrationNamespace(scenario);

  std::vector<std::string> arguments;
  arguments.reserve(6 + parameter_overrides.size() * 2);
  arguments.emplace_back(RAYSTAR_NODE_PATH);
  arguments.emplace_back("--ros-args");
  arguments.emplace_back("-r");
  arguments.emplace_back("__ns:=" + spawned.node_namespace);
  // map_topic is read-only in the node.  Remapping it here keeps the global
  // /map name out of the test domain even when another map publisher is live.
  arguments.emplace_back("-p");
  arguments.emplace_back("map_topic:=" + spawned.node_namespace + "/map");
  for (const auto& override : parameter_overrides)
  {
    arguments.emplace_back("-p");
    arguments.push_back(override);
  }

  int error_pipe[2] = {-1, -1};
  if (pipe(error_pipe) != 0)
  {
    spawned.launch_errno = errno;
    return spawned;
  }
  const int flags = fcntl(error_pipe[1], F_GETFD);
  if (flags < 0 || fcntl(error_pipe[1], F_SETFD, flags | FD_CLOEXEC) < 0)
  {
    spawned.launch_errno = errno;
    close(error_pipe[0]);
    close(error_pipe[1]);
    return spawned;
  }

  std::vector<char *> argv;
  argv.reserve(arguments.size() + 1);
  for (auto& argument : arguments)
    argv.push_back(argument.data());
  argv.push_back(nullptr);

  spawned.pid = fork();
  if (spawned.pid < 0)
  {
    spawned.launch_errno = errno;
    close(error_pipe[0]);
    close(error_pipe[1]);
    return spawned;
  }
  if (spawned.pid == 0)
  {
    close(error_pipe[0]);
    execv(argv[0], argv.data());
    const int child_errno = errno;
    // write(2) is async-signal-safe and therefore suitable in the post-fork
    // child even though the parent process owns DDS/background threads.
    const char *bytes = reinterpret_cast<const char *>(&child_errno);
    size_t written = 0;
    while (written < sizeof(child_errno))
    {
      const ssize_t count = write(
        error_pipe[1], bytes + written, sizeof(child_errno) - written);
      if (count > 0)
        written += static_cast<size_t>(count);
      else if (count < 0 && errno == EINTR)
        continue;
      else
        break;
    }
    _exit(127);
  }

  close(error_pipe[1]);
  int child_errno = 0;
  size_t received = 0;
  while (received < sizeof(child_errno))
  {
    const ssize_t count = read(
      error_pipe[0], reinterpret_cast<char *>(&child_errno) + received,
      sizeof(child_errno) - received);
    if (count > 0)
      received += static_cast<size_t>(count);
    else if (count < 0 && errno == EINTR)
      continue;
    else
      break;
  }
  close(error_pipe[0]);
  if (received != 0)
  {
    spawned.launch_errno = child_errno;
    int ignored_status = 0;
    (void)waitpid(spawned.pid, &ignored_status, 0);
    spawned.pid = -1;
  }
  return spawned;
}

static bool waitForChildExit(
  pid_t pid, int& status, std::chrono::milliseconds timeout)
{
  const auto deadline = std::chrono::steady_clock::now() + timeout;
  while (std::chrono::steady_clock::now() < deadline)
  {
    const pid_t result = waitpid(pid, &status, WNOHANG);
    if (result == pid)
      return true;
    if (result < 0 && errno != EINTR)
      return false;
    std::this_thread::sleep_for(10ms);
  }
  return waitpid(pid, &status, WNOHANG) == pid;
}

static RaystarService::Response::SharedPtr callService(
  rclcpp::executors::SingleThreadedExecutor& executor,
  const rclcpp::Client<RaystarService>::SharedPtr& client,
  const RaystarService::Request::SharedPtr& request)
{
  auto future = client->async_send_request(request);
  const auto end_time = std::chrono::steady_clock::now() + 10s;
  while (std::chrono::steady_clock::now() < end_time)
  {
    if (future.wait_for(100ms) == std::future_status::ready)
      break;
    executor.spin_some(100ms);
  }
  if (future.wait_for(0s) != std::future_status::ready)
  {
    ADD_FAILURE() << "Service call did not complete within 10 seconds";
    return nullptr;
  }
  return future.get();
}

static std::string lowercase(std::string text)
{
  std::transform(text.begin(), text.end(), text.begin(),
    [](unsigned char character) {
      return static_cast<char>(std::tolower(character));
    });
  return text;
}

template<typename ResponseT>
static void expectStructuredResultInvariants(const ResponseT& response)
{
  const auto& info = response.result_info;
  EXPECT_EQ(info.returned_path_count, response.path_results.size());
  EXPECT_LE(info.returned_path_count, info.found_path_count);
  EXPECT_LE(info.found_path_count, info.requested_path_count);
  EXPECT_EQ(response.success, !response.path_results.empty());
  if (info.request_satisfied)
  {
    EXPECT_TRUE(info.search_complete);
    EXPECT_EQ(info.returned_path_count, info.requested_path_count);
  }
}

static void expectIndependentCollisionFreePath(
  const nav_msgs::msg::OccupancyGrid& map, const nav_msgs::msg::Path& path,
  bool allow_unknown, int occupied_threshold = 99)
{
  raystar::test_oracle::OracleOptions options;
  options.occupied_threshold = occupied_threshold;
  options.allow_unknown = allow_unknown;
  const auto result = raystar::test_oracle::validatePath(map, path, options);
  EXPECT_EQ(result.status,
    raystar::test_oracle::ValidationStatus::kCollisionFree)
    << result.diagnostic;
}

static void expectIndependentSelfIntersectionFreePath(
  const nav_msgs::msg::Path& path)
{
  const auto result =
    raystar::test_oracle::validateNoSelfIntersection(path);
  EXPECT_EQ(result.status,
    raystar::test_oracle::SelfIntersectionStatus::kIntersectionFree)
    << result.diagnostic;
}

static RaystarService::Request::SharedPtr makeSmallGridPropertyRequest(
  const raystar::test_property::SmallGridCase& property_case)
{
  auto request = std::make_shared<RaystarService::Request>();
  request->map = property_case.map;
  request->start.header.frame_id = property_case.map.header.frame_id;
  request->goal.header.frame_id = property_case.map.header.frame_id;

  const double resolution =
    static_cast<double>(property_case.map.info.resolution);
  const auto cell_center = [&](const raystar::test_property::GridCell& cell) {
    return std::pair<double, double>{
      property_case.map.info.origin.position.x +
        (static_cast<double>(cell.x) + 0.5) * resolution,
      property_case.map.info.origin.position.y +
        (static_cast<double>(cell.y) + 0.5) * resolution};
  };
  const auto start = cell_center(property_case.start);
  const auto goal = cell_center(property_case.goal);
  request->start.pose.position.x = start.first;
  request->start.pose.position.y = start.second;
  request->start.pose.orientation.w = 1.0;
  request->goal.pose.position.x = goal.first;
  request->goal.pose.position.y = goal.second;
  request->goal.pose.orientation.w = 1.0;
  request->k = 1;
  request->allow_self_crossing = false;
  request->allow_unknown = property_case.allow_unknown;
  request->include_debug = false;
  return request;
}

static double polylineMetricLength(const nav_msgs::msg::Path& path)
{
  double length = 0.0;
  for (size_t index = 1; index < path.poses.size(); ++index)
  {
    const auto& previous = path.poses[index - 1].pose.position;
    const auto& current = path.poses[index].pose.position;
    length += std::hypot(
      current.x - previous.x, current.y - previous.y);
  }
  return length;
}

static void expectSmallGridPathContract(
  const raystar::test_property::SmallGridCase& property_case,
  const RaystarService::Request& request,
  const raystar_interfaces::msg::PathResult& path_result)
{
  const auto& path = path_result.path;
  EXPECT_EQ(path.header.frame_id, property_case.map.header.frame_id);
  ASSERT_GE(path.poses.size(), 2u);
  for (size_t pose_index = 0; pose_index < path.poses.size(); ++pose_index)
  {
    SCOPED_TRACE("pose index " + std::to_string(pose_index));
    const auto& pose_stamped = path.poses[pose_index];
    const auto& position = pose_stamped.pose.position;
    const auto& orientation = pose_stamped.pose.orientation;
    EXPECT_EQ(pose_stamped.header.frame_id,
      property_case.map.header.frame_id);
    EXPECT_TRUE(std::isfinite(position.x));
    EXPECT_TRUE(std::isfinite(position.y));
    EXPECT_TRUE(std::isfinite(position.z));
    EXPECT_DOUBLE_EQ(position.z, 0.0);
    EXPECT_TRUE(std::isfinite(orientation.x));
    EXPECT_TRUE(std::isfinite(orientation.y));
    EXPECT_TRUE(std::isfinite(orientation.z));
    EXPECT_TRUE(std::isfinite(orientation.w));
    EXPECT_DOUBLE_EQ(orientation.x, 0.0);
    EXPECT_DOUBLE_EQ(orientation.y, 0.0);
    EXPECT_DOUBLE_EQ(orientation.z, 0.0);
    EXPECT_DOUBLE_EQ(orientation.w, 1.0);
  }

  const auto& first = path.poses.front().pose.position;
  const auto& last = path.poses.back().pose.position;
  EXPECT_NEAR(first.x, request.start.pose.position.x, 1e-9);
  EXPECT_NEAR(first.y, request.start.pose.position.y, 1e-9);
  EXPECT_NEAR(last.x, request.goal.pose.position.x, 1e-9);
  EXPECT_NEAR(last.y, request.goal.pose.position.y, 1e-9);

  EXPECT_TRUE(std::isfinite(path_result.cost));
  EXPECT_GE(path_result.cost, 0.0);
  const double sampled_length = polylineMetricLength(path);
  const double length_tolerance = 1e-8 * std::max(1.0, sampled_length);
  EXPECT_NEAR(path_result.cost, sampled_length, length_tolerance);
  const double direct_distance = std::hypot(
    request.goal.pose.position.x - request.start.pose.position.x,
    request.goal.pose.position.y - request.start.pose.position.y);
  EXPECT_GE(path_result.cost + length_tolerance, direct_distance);
  expectIndependentCollisionFreePath(
    property_case.map, path, property_case.allow_unknown,
    raystar::test_property::kOccupiedThreshold);
  expectIndependentSelfIntersectionFreePath(path);
}

static void expectSmallGridReachableResult(
  const raystar::test_property::SmallGridCase& property_case,
  const RaystarService::Request& request,
  const RaystarService::Response& response)
{
  const auto& info = response.result_info;
  EXPECT_TRUE(response.success) << response.message;
  EXPECT_EQ(info.status, PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_EQ(info.limits_reached, PlanningResultInfo::LIMIT_NONE);
  EXPECT_TRUE(info.request_satisfied);
  EXPECT_TRUE(info.search_complete);
  EXPECT_TRUE(info.output_complete);
  EXPECT_FALSE(info.debug_requested);
  EXPECT_TRUE(info.debug_output_complete);
  EXPECT_EQ(info.requested_path_count, 1u);
  EXPECT_EQ(info.found_path_count, 1u);
  EXPECT_EQ(info.returned_path_count, 1u);
  EXPECT_TRUE(raystar_interfaces::mapIdsEqual(
    info.map_id, raystar_interfaces::computeMapId(property_case.map)));
  EXPECT_TRUE(std::isfinite(info.map_time_ms));
  EXPECT_GE(info.map_time_ms, 0.0);
  EXPECT_TRUE(std::isfinite(info.plan_time_ms));
  EXPECT_GE(info.plan_time_ms, 0.0);
  EXPECT_TRUE(response.debug_nodes.empty());
  ASSERT_EQ(response.path_results.size(), 1u);
  expectSmallGridPathContract(
    property_case, request, response.path_results.front());
  expectStructuredResultInvariants(response);
}

static void expectSmallGridUnreachableResult(
  const raystar::test_property::SmallGridCase& property_case,
  const RaystarService::Response& response)
{
  const auto& info = response.result_info;
  EXPECT_FALSE(response.success);
  EXPECT_EQ(info.status, PlanningResultInfo::STATUS_NO_PATH);
  EXPECT_EQ(info.limits_reached, PlanningResultInfo::LIMIT_NONE);
  EXPECT_FALSE(info.request_satisfied);
  EXPECT_TRUE(info.search_complete);
  EXPECT_TRUE(info.output_complete);
  EXPECT_FALSE(info.debug_requested);
  EXPECT_TRUE(info.debug_output_complete);
  EXPECT_EQ(info.requested_path_count, 1u);
  EXPECT_EQ(info.found_path_count, 0u);
  EXPECT_EQ(info.returned_path_count, 0u);
  EXPECT_TRUE(raystar_interfaces::mapIdsEqual(
    info.map_id, raystar_interfaces::computeMapId(property_case.map)));
  EXPECT_TRUE(std::isfinite(info.map_time_ms));
  EXPECT_GE(info.map_time_ms, 0.0);
  EXPECT_TRUE(std::isfinite(info.plan_time_ms));
  EXPECT_GE(info.plan_time_ms, 0.0);
  EXPECT_TRUE(response.path_results.empty());
  EXPECT_TRUE(response.debug_nodes.empty());
  expectStructuredResultInvariants(response);
}

static void expectStableSmallGridResponses(
  const RaystarService::Response& first,
  const RaystarService::Response& second)
{
  EXPECT_EQ(first.success, second.success);
  const auto& first_info = first.result_info;
  const auto& second_info = second.result_info;
  EXPECT_TRUE(raystar_interfaces::mapIdsEqual(
    first_info.map_id, second_info.map_id));
  EXPECT_EQ(first_info.status, second_info.status);
  EXPECT_EQ(first_info.limits_reached, second_info.limits_reached);
  EXPECT_EQ(first_info.request_satisfied, second_info.request_satisfied);
  EXPECT_EQ(first_info.search_complete, second_info.search_complete);
  EXPECT_EQ(first_info.output_complete, second_info.output_complete);
  EXPECT_EQ(first_info.debug_requested, second_info.debug_requested);
  EXPECT_EQ(first_info.debug_output_complete,
    second_info.debug_output_complete);
  EXPECT_EQ(first_info.requested_path_count,
    second_info.requested_path_count);
  EXPECT_EQ(first_info.found_path_count, second_info.found_path_count);
  EXPECT_EQ(first_info.returned_path_count,
    second_info.returned_path_count);
  EXPECT_EQ(first.path_results.size(), second.path_results.size());
  EXPECT_EQ(first.debug_nodes.size(), second.debug_nodes.size());

  ASSERT_EQ(first.path_results.size(), second.path_results.size());
  for (size_t path_index = 0;
       path_index < first.path_results.size(); ++path_index)
  {
    SCOPED_TRACE("repeated path index " + std::to_string(path_index));
    const auto& first_result = first.path_results[path_index];
    const auto& second_result = second.path_results[path_index];
    EXPECT_DOUBLE_EQ(first_result.cost, second_result.cost);
    EXPECT_EQ(first_result.path.header.frame_id,
      second_result.path.header.frame_id);
    ASSERT_EQ(first_result.path.poses.size(),
      second_result.path.poses.size());
    for (size_t pose_index = 0;
         pose_index < first_result.path.poses.size(); ++pose_index)
    {
      SCOPED_TRACE("repeated pose index " + std::to_string(pose_index));
      const auto& first_pose = first_result.path.poses[pose_index];
      const auto& second_pose = second_result.path.poses[pose_index];
      EXPECT_EQ(first_pose.header.frame_id, second_pose.header.frame_id);
      EXPECT_DOUBLE_EQ(first_pose.pose.position.x,
        second_pose.pose.position.x);
      EXPECT_DOUBLE_EQ(first_pose.pose.position.y,
        second_pose.pose.position.y);
      EXPECT_DOUBLE_EQ(first_pose.pose.position.z,
        second_pose.pose.position.z);
      EXPECT_DOUBLE_EQ(first_pose.pose.orientation.x,
        second_pose.pose.orientation.x);
      EXPECT_DOUBLE_EQ(first_pose.pose.orientation.y,
        second_pose.pose.orientation.y);
      EXPECT_DOUBLE_EQ(first_pose.pose.orientation.z,
        second_pose.pose.orientation.z);
      EXPECT_DOUBLE_EQ(first_pose.pose.orientation.w,
        second_pose.pose.orientation.w);
    }
  }
}

template<typename ResponseT>
static void expectSingleObstacleSimplePathContract(
  const nav_msgs::msg::OccupancyGrid& map, const ResponseT& response)
{
  const auto& info = response.result_info;
  EXPECT_TRUE(response.success) << response.message;
  EXPECT_EQ(info.status, PlanningResultInfo::STATUS_FEWER_PATHS);
  EXPECT_EQ(info.limits_reached, PlanningResultInfo::LIMIT_NONE);
  EXPECT_EQ(info.requested_path_count, 3u);
  EXPECT_EQ(info.found_path_count, 2u);
  EXPECT_EQ(info.returned_path_count, 2u);
  EXPECT_TRUE(info.search_complete);
  EXPECT_TRUE(info.output_complete);
  EXPECT_FALSE(info.request_satisfied);
  EXPECT_FALSE(info.debug_requested);
  EXPECT_TRUE(info.debug_output_complete);
  EXPECT_TRUE(raystar_interfaces::mapIdsEqual(
    info.map_id, raystar_interfaces::computeMapId(map)));
  EXPECT_TRUE(std::isfinite(info.map_time_ms));
  EXPECT_GE(info.map_time_ms, 0.0);
  EXPECT_TRUE(std::isfinite(info.plan_time_ms));
  EXPECT_GE(info.plan_time_ms, 0.0);
  EXPECT_GT(info.expanded_nodes, 0u);
  EXPECT_EQ(response.path_results.size(), 2u);
  expectStructuredResultInvariants(response);

  for (size_t index = 0; index < response.path_results.size(); ++index)
  {
    SCOPED_TRACE("path index " + std::to_string(index));
    const auto& path = response.path_results[index].path;
    expectIndependentCollisionFreePath(map, path, false);
    expectIndependentSelfIntersectionFreePath(path);
  }
}

class ChildProcessGuard
{
public:
  explicit ChildProcessGuard(pid_t pid) : pid_(pid) {}
  ChildProcessGuard(const ChildProcessGuard&) = delete;
  ChildProcessGuard& operator=(const ChildProcessGuard&) = delete;

  ~ChildProcessGuard()
  {
    stopAndValidate();
  }

  void disarm() noexcept
  {
    pid_ = -1;
  }

private:
  void stopAndValidate() noexcept
  {
    if (pid_ <= 0)
      return;

    int status = 0;
    bool reaped = waitForChildExit(pid_, status, 0ms);
    bool used_sigterm = false;
    bool used_sigkill = false;
    if (!reaped)
    {
      if (kill(pid_, SIGINT) < 0 && errno != ESRCH)
        ADD_FAILURE() << "Failed to send SIGINT to raystar_node: "
                      << std::strerror(errno);
      reaped = waitForChildExit(pid_, status, 3s);
    }
    if (!reaped)
    {
      used_sigterm = true;
      if (kill(pid_, SIGTERM) < 0 && errno != ESRCH)
        ADD_FAILURE() << "Failed to send SIGTERM to raystar_node: "
                      << std::strerror(errno);
      reaped = waitForChildExit(pid_, status, 2s);
    }
    if (!reaped)
    {
      used_sigkill = true;
      if (kill(pid_, SIGKILL) < 0 && errno != ESRCH)
        ADD_FAILURE() << "Failed to send SIGKILL to raystar_node: "
                      << std::strerror(errno);
      reaped = waitForChildExit(pid_, status, 2s);
    }

    if (!reaped)
    {
      ADD_FAILURE() << "raystar_node child did not exit within bounded "
                    << "SIGINT/SIGTERM/SIGKILL cleanup";
      return;
    }
    if (used_sigkill)
      ADD_FAILURE() << "raystar_node required SIGKILL during cleanup";
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0)
    {
      ADD_FAILURE() << "raystar_node exited unexpectedly (status=" << status
                    << ")";
    }
    pid_ = -1;
    (void)used_sigterm;
  }

  pid_t pid_;
};

class IntegrationTestFixture : public ::testing::Test
{
protected:
  static void SetUpTestSuite()
  {
    // Keep DDS discovery and rosout files isolated from a concurrently
    // running integration binary.  A caller-provided domain is retained so
    // CI can deliberately place tests in a preselected domain; otherwise a
    // PID-derived valid domain is assigned before rclcpp initializes.
    const auto existing_domain = std::getenv("ROS_DOMAIN_ID");
    const int domain_id = existing_domain ? std::atoi(existing_domain) :
      20 + static_cast<int>(getpid() % 180);
    ASSERT_GE(domain_id, 0);
    ASSERT_LE(domain_id, 232);
    const std::string domain_text = std::to_string(domain_id);
    ASSERT_EQ(setenv("ROS_DOMAIN_ID", domain_text.c_str(), 1), 0)
      << std::strerror(errno);

    const std::filesystem::path log_directory =
      std::filesystem::temp_directory_path() /
      ("raystar_integration_" + std::to_string(static_cast<long long>(getpid())));
    std::error_code directory_error;
    std::filesystem::create_directories(log_directory, directory_error);
    ASSERT_FALSE(directory_error)
      << "Unable to create ROS_LOG_DIR " << log_directory << ": "
      << directory_error.message();
    ASSERT_EQ(setenv("ROS_LOG_DIR", log_directory.c_str(), 1), 0)
      << std::strerror(errno);

    rclcpp::init(0, nullptr);

    const auto spawned = spawnRaystarNode("main");
    ASSERT_GT(spawned.pid, 0)
      << "Unable to exec raystar_node: " <<
      (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
      "fork/pipe failure");
    pid_ = spawned.pid;
    main_namespace_ = spawned.node_namespace;

    auto node = std::make_shared<rclcpp::Node>("test_wait_node");
    auto client = node->create_client<raystar_interfaces::srv::GetRaystarPaths>(
      raystarEndpoint(main_namespace_, "get_raystar_paths"));
    auto action_client = rclcpp_action::create_client<RaystarAction>(
      node, raystarEndpoint(main_namespace_, "plan_paths"));

    bool found = false;
    for (int i = 0; i < 50; ++i) {
      if (client->wait_for_service(100ms)) { found = true; break; }
    }
    ASSERT_TRUE(found) << "Service "
      << raystarEndpoint(main_namespace_, "get_raystar_paths")
      << " not available after 5s";
    ASSERT_TRUE(action_client->wait_for_action_server(5s))
      << "Action " << raystarEndpoint(main_namespace_, "plan_paths")
      << " not available after 5s";

    rclcpp::executors::SingleThreadedExecutor executor;
    executor.add_node(node);
    ASSERT_TRUE(cacheMapAndWait(executor, node, makeTestGrid(), main_namespace_))
      << "Raystar server did not admit the default cached map";
  }

  static void TearDownTestSuite()
  {
    if (pid_ > 0)
    {
      ChildProcessGuard process(pid_);
      pid_ = -1;
    }
    auto ctx = rclcpp::contexts::get_global_default_context();
    if (ctx && ctx->is_valid()) {
      ctx->shutdown("test teardown");
    }
  }

  static pid_t pid_;
  static std::string main_namespace_;
};

pid_t IntegrationTestFixture::pid_ = 0;
std::string IntegrationTestFixture::main_namespace_ = mainServerNamespace();

TEST_F(IntegrationTestFixture, ServiceCallReturnsPaths)
{
  auto node = std::make_shared<rclcpp::Node>("test_integration_client");
  auto client = node->create_client<raystar_interfaces::srv::GetRaystarPaths>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));

  ASSERT_TRUE(client->wait_for_service(2s));

  auto request = makeTestRequest();
  request->k = 1;

  auto result_future = client->async_send_request(request);

  // Spin with a single-threaded executor until the result is ready
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);
  auto end_time = std::chrono::steady_clock::now() + 10s;
  while (std::chrono::steady_clock::now() < end_time) {
    auto status = result_future.wait_for(100ms);
    if (status == std::future_status::ready) break;
    executor.spin_some(100ms);
  }
  ASSERT_EQ(result_future.wait_for(0s), std::future_status::ready);

  auto result = result_future.get();
  EXPECT_TRUE(result->success);
  EXPECT_GE(result->path_results.size(), 1u);
  EXPECT_EQ(result->result_info.status, PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_TRUE(result->result_info.request_satisfied);
  EXPECT_TRUE(result->result_info.search_complete);
  EXPECT_TRUE(result->result_info.output_complete);
  expectStructuredResultInvariants(*result);
}

TEST_F(IntegrationTestFixture, ExhaustedSearchReportsFewerPathsWithoutStringParsing)
{
  auto node = std::make_shared<rclcpp::Node>("test_fewer_paths_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto request = makeTestRequest();
  std::fill(request->map.data.begin(), request->map.data.end(), 0);
  request->k = 3;
  const auto response = callService(executor, client, request);

  ASSERT_NE(response, nullptr);
  ASSERT_TRUE(response->success) << response->message;
  ASSERT_EQ(response->path_results.size(), 1u);
  EXPECT_EQ(response->result_info.status,
    PlanningResultInfo::STATUS_FEWER_PATHS);
  EXPECT_EQ(response->result_info.requested_path_count, 3u);
  EXPECT_EQ(response->result_info.found_path_count, 1u);
  EXPECT_EQ(response->result_info.returned_path_count, 1u);
  EXPECT_TRUE(response->result_info.search_complete);
  EXPECT_TRUE(response->result_info.output_complete);
  EXPECT_FALSE(response->result_info.request_satisfied);
  expectStructuredResultInvariants(*response);
}

TEST_F(IntegrationTestFixture,
  SelfCrossingPolicyIsEnforcedAcrossServiceAndActionOutputs)
{
  auto node = std::make_shared<rclcpp::Node>(
    "test_self_crossing_contract_client");
  auto service_client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));
  auto action_client = rclcpp_action::create_client<RaystarAction>(
    node, raystarEndpoint(main_namespace_, "plan_paths"));
  ASSERT_TRUE(service_client->wait_for_service(2s));
  ASSERT_TRUE(action_client->wait_for_action_server(5s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  // A single rectangular obstacle has exactly two simple path classes for
  // these endpoints.  With K=3 and self-crossing disabled, normal frontier
  // exhaustion must therefore return the upper and lower simple paths as a
  // complete, untruncated FEWER_PATHS result.
  auto service_request = makeTestRequest();
  service_request->k = 3;
  service_request->allow_self_crossing = false;
  const auto service_response = callService(
    executor, service_client, service_request);
  ASSERT_NE(service_response, nullptr);
  expectSingleObstacleSimplePathContract(
    service_request->map, *service_response);

  ASSERT_TRUE(cacheMapAndWait(
      executor, node, service_request->map, main_namespace_))
    << "Raystar server did not admit the self-crossing contract map";
  auto action_goal = makeTestActionGoal();
  action_goal.k = 3;
  action_goal.allow_self_crossing = false;
  auto goal_future = action_client->async_send_goal(action_goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s))
    << "Self-crossing-disabled Action goal was not acknowledged";
  const RaystarGoalHandle::SharedPtr goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr)
    << "Self-crossing-disabled Action goal was rejected";
  auto result_future = action_client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s))
    << "Self-crossing-disabled Action goal did not complete";
  const auto action_result = result_future.get();
  EXPECT_EQ(action_result.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(action_result.result, nullptr);
  expectSingleObstacleSimplePathContract(
    service_request->map, *action_result.result);

  // Enabling the policy only permits self-crossing candidates; it does not
  // require the planner to manufacture one.  This smoke request therefore
  // verifies acceptance and the ordinary output contract without asserting a
  // particular self-intersecting path or a policy-dependent path count.
  auto permissive_goal = makeTestActionGoal();
  permissive_goal.k = 1;
  permissive_goal.allow_self_crossing = true;
  auto permissive_goal_future =
    action_client->async_send_goal(permissive_goal);
  ASSERT_TRUE(waitForFuture(executor, permissive_goal_future, 5s))
    << "Self-crossing-enabled Action goal was not acknowledged";
  const RaystarGoalHandle::SharedPtr permissive_handle =
    permissive_goal_future.get();
  ASSERT_NE(permissive_handle, nullptr)
    << "Self-crossing-enabled Action goal was rejected";
  auto permissive_result_future =
    action_client->async_get_result(permissive_handle);
  ASSERT_TRUE(waitForFuture(executor, permissive_result_future, 10s))
    << "Self-crossing-enabled Action goal did not complete";
  const auto permissive_result = permissive_result_future.get();
  EXPECT_EQ(permissive_result.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(permissive_result.result, nullptr);
  const auto& permissive_response = *permissive_result.result;
  EXPECT_TRUE(permissive_response.success) << permissive_response.message;
  EXPECT_EQ(permissive_response.result_info.status,
    PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_EQ(permissive_response.result_info.limits_reached,
    PlanningResultInfo::LIMIT_NONE);
  EXPECT_EQ(permissive_response.result_info.requested_path_count, 1u);
  EXPECT_EQ(permissive_response.result_info.found_path_count, 1u);
  EXPECT_EQ(permissive_response.result_info.returned_path_count, 1u);
  EXPECT_TRUE(permissive_response.result_info.request_satisfied);
  EXPECT_TRUE(permissive_response.result_info.search_complete);
  EXPECT_TRUE(permissive_response.result_info.output_complete);
  EXPECT_FALSE(permissive_response.result_info.debug_requested);
  EXPECT_TRUE(permissive_response.result_info.debug_output_complete);
  EXPECT_TRUE(raystar_interfaces::mapIdsEqual(
    permissive_response.result_info.map_id,
    raystar_interfaces::computeMapId(service_request->map)));
  ASSERT_EQ(permissive_response.path_results.size(), 1u);
  expectIndependentCollisionFreePath(
    service_request->map,
    permissive_response.path_results.front().path, false);
  expectStructuredResultInvariants(permissive_response);
}

TEST_F(IntegrationTestFixture, OccupancyThresholdAndUnknownPolicyAreApplied)
{
  auto node = std::make_shared<rclcpp::Node>(
    "test_occupancy_threshold_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto below_threshold = callService(
    executor, client, makeVerticalBarrierRequest(98));
  ASSERT_NE(below_threshold, nullptr);
  EXPECT_TRUE(below_threshold->success) << below_threshold->message;

  auto at_threshold = callService(
    executor, client, makeVerticalBarrierRequest(99));
  ASSERT_NE(at_threshold, nullptr);
  EXPECT_FALSE(at_threshold->success);
  EXPECT_NE(at_threshold->message.find("No path exists"), std::string::npos)
    << at_threshold->message;

  auto unknown_blocked = callService(
    executor, client, makeVerticalBarrierRequest(-1, false));
  ASSERT_NE(unknown_blocked, nullptr);
  EXPECT_FALSE(unknown_blocked->success);
  EXPECT_NE(unknown_blocked->message.find("No path exists"), std::string::npos)
    << unknown_blocked->message;

  auto unknown_allowed = callService(
    executor, client, makeVerticalBarrierRequest(-1, true));
  ASSERT_NE(unknown_allowed, nullptr);
  EXPECT_TRUE(unknown_allowed->success) << unknown_allowed->message;
}

TEST_F(IntegrationTestFixture, FixedSeedSmallGridDifferentialProperties)
{
  auto node = std::make_shared<rclcpp::Node>(
    "test_small_grid_property_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  constexpr size_t property_case_count = 24;
  for (size_t case_index = 0;
       case_index < property_case_count; ++case_index)
  {
    const auto property_case =
      raystar::test_property::makeSmallGridCase(case_index);
    SCOPED_TRACE(raystar::test_property::describeSmallGridCase(property_case));

    const bool bfs_reachable =
      raystar::test_property::fourNeighborReachable(property_case);
    ASSERT_EQ(bfs_reachable, property_case.expected_reachable)
      << "The fixed-seed generator did not preserve its reachability label";
    ASSERT_FALSE(
      raystar::test_property::hasUnsupportedDiagonalContact(property_case))
      << "The fixed-seed generator produced unsupported contour topology";

    const auto request = makeSmallGridPropertyRequest(property_case);
    const auto first_response = callService(executor, client, request);
    ASSERT_NE(first_response, nullptr);
    const auto second_response = callService(executor, client, request);
    ASSERT_NE(second_response, nullptr);

    EXPECT_EQ(first_response->success, bfs_reachable)
      << "Raystar reachability differs from the independent four-neighbor BFS";
    EXPECT_EQ(second_response->success, bfs_reachable)
      << "Repeated Raystar reachability differs from the independent BFS";
    if (bfs_reachable)
    {
      expectSmallGridReachableResult(
        property_case, *request, *first_response);
      expectSmallGridReachableResult(
        property_case, *request, *second_response);
    }
    else
    {
      expectSmallGridUnreachableResult(property_case, *first_response);
      expectSmallGridUnreachableResult(property_case, *second_response);
    }
    expectStableSmallGridResponses(*first_response, *second_response);
  }
}

TEST_F(IntegrationTestFixture, ActionCancelRejectsConcurrentGoalAndRecovers)
{
  auto node = std::make_shared<rclcpp::Node>("test_action_cancel_client");
  auto client = rclcpp_action::create_client<RaystarAction>(
    node, raystarEndpoint(main_namespace_, "plan_paths"));
  auto service_client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_action_server(5s));
  ASSERT_TRUE(service_client->wait_for_service(5s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  ASSERT_TRUE(cacheMapAndWait(
      executor, node, makeLongRunningActionMap(), main_namespace_))
    << "The long-running map was not admitted by the Action server";
  auto first_goal_future = client->async_send_goal(makeLongRunningActionGoal());
  ASSERT_TRUE(waitForFuture(executor, first_goal_future, 10s))
    << "The long-running Action goal was not acknowledged";
  const RaystarGoalHandle::SharedPtr first_handle = first_goal_future.get();
  ASSERT_NE(first_handle, nullptr) << "The first Action goal was rejected";
  auto first_result_future = client->async_get_result(first_handle);

  // Goal admission is serialized by the server.  Waiting for the first goal
  // handle before sending this request makes the assertion independent of
  // discovery/transport ordering: the first goal is already the active goal.
  auto second_goal_future = client->async_send_goal(makeTestActionGoal());
  const bool second_response_ready =
    waitForFuture(executor, second_goal_future, 5s);
  RaystarGoalHandle::SharedPtr second_handle;
  if (second_response_ready)
    second_handle = second_goal_future.get();
  else
    ADD_FAILURE() << "The concurrent Action goal received no admission response";

  auto busy_service_future = service_client->async_send_request(makeTestRequest());
  const bool busy_service_ready =
    waitForFuture(executor, busy_service_future, 5s);
  RaystarService::Response::SharedPtr busy_service_response;
  if (busy_service_ready)
    busy_service_response = busy_service_future.get();
  else
    ADD_FAILURE() << "The compatibility Service received no busy response";

  auto cancel_future = client->async_cancel_goal(first_handle);
  const bool cancel_response_ready = waitForFuture(executor, cancel_future, 5s);
  if (!cancel_response_ready)
  {
    ADD_FAILURE() << "Cancel request received no response";
  }
  else
  {
    const auto cancel_response = cancel_future.get();
    EXPECT_FALSE(cancel_response->goals_canceling.empty())
      << "The active goal was not accepted for cancellation";
  }

  // If a broken server accepted the second goal, request its cancellation as
  // cleanup before making the rejection assertion below.
  if (second_handle)
    (void)client->async_cancel_goal(second_handle);

  const bool canceled_result_ready =
    waitForFuture(executor, first_result_future, 15s);
  ASSERT_TRUE(canceled_result_ready)
    << "Canceled Action goal did not reach a terminal state";
  const auto canceled_result = first_result_future.get();

  EXPECT_TRUE(second_response_ready);
  EXPECT_EQ(second_handle, nullptr)
    << "A second Action goal was accepted while planning was active";
  EXPECT_TRUE(busy_service_ready);
  ASSERT_NE(busy_service_response, nullptr);
  EXPECT_FALSE(busy_service_response->success);
  EXPECT_EQ(busy_service_response->result_info.status,
    PlanningResultInfo::STATUS_BUSY);
  EXPECT_EQ(busy_service_response->result_info.requested_path_count, 3u);
  expectStructuredResultInvariants(*busy_service_response);
  EXPECT_NE(lowercase(busy_service_response->message).find("busy"),
    std::string::npos) << busy_service_response->message;
  EXPECT_EQ(canceled_result.code, rclcpp_action::ResultCode::CANCELED);
  ASSERT_NE(canceled_result.result, nullptr);
  EXPECT_FALSE(canceled_result.result->success);
  EXPECT_EQ(canceled_result.result->result_info.status,
    PlanningResultInfo::STATUS_CANCELLED);
  EXPECT_NE(canceled_result.result->result_info.limits_reached &
    PlanningResultInfo::LIMIT_CANCELLED, 0u);
  EXPECT_FALSE(canceled_result.result->result_info.search_complete);
  EXPECT_FALSE(canceled_result.result->result_info.output_complete);
  EXPECT_FALSE(canceled_result.result->result_info.request_satisfied);
  EXPECT_TRUE(canceled_result.result->path_results.empty());
  expectStructuredResultInvariants(*canceled_result.result);
  EXPECT_NE(lowercase(canceled_result.result->message).find("cancel"),
    std::string::npos) << canceled_result.result->message;

  ASSERT_TRUE(cacheMapAndWait(
      executor, node, makeTestGrid(), main_namespace_))
    << "The recovery map was not admitted by the Action server";
  auto recovery_goal_future = client->async_send_goal(makeTestActionGoal());
  ASSERT_TRUE(waitForFuture(executor, recovery_goal_future, 5s))
    << "The recovery Action goal was not acknowledged";
  const RaystarGoalHandle::SharedPtr recovery_handle = recovery_goal_future.get();
  ASSERT_NE(recovery_handle, nullptr)
    << "The server did not admit a goal after cancellation completed";

  auto recovery_result_future = client->async_get_result(recovery_handle);
  ASSERT_TRUE(waitForFuture(executor, recovery_result_future, 10s))
    << "The recovery Action goal did not complete";
  const auto recovery_result = recovery_result_future.get();
  EXPECT_EQ(recovery_result.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(recovery_result.result, nullptr);
  EXPECT_TRUE(recovery_result.result->success)
    << recovery_result.result->message;
  EXPECT_FALSE(recovery_result.result->path_results.empty());
  EXPECT_TRUE(recovery_result.result->result_info.search_complete);
  EXPECT_TRUE(recovery_result.result->result_info.output_complete);
  expectStructuredResultInvariants(*recovery_result.result);
}

TEST_F(IntegrationTestFixture, SigintDuringActiveActionJoinsWorkerAndExits)
{
  const auto spawned = spawnRaystarNode(
    "graceful_shutdown", {"planning_timeout_ms:=60000"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec graceful-shutdown raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  const pid_t shutdown_pid = spawned.pid;
  ChildProcessGuard shutdown_process(shutdown_pid);

  auto node = std::make_shared<rclcpp::Node>("test_action_shutdown_client");
  auto client = rclcpp_action::create_client<RaystarAction>(
    node, raystarEndpoint(spawned.node_namespace, "plan_paths"));
  ASSERT_TRUE(client->wait_for_action_server(5s));

  bool worker_started = false;
  std::size_t deleteall_count = 0;
  auto worker_subscription =
    node->create_subscription<visualization_msgs::msg::MarkerArray>(
      raystarEndpoint(spawned.node_namespace, "non_homotopic_paths"),
      rclcpp::QoS(rclcpp::KeepLast(10)).transient_local(),
      [&](const visualization_msgs::msg::MarkerArray::ConstSharedPtr message) {
        const auto count = static_cast<std::size_t>(std::count_if(
          message->markers.begin(), message->markers.end(),
          [](const auto& marker) {
            return marker.action == visualization_msgs::msg::Marker::DELETEALL;
          }));
        deleteall_count += count;
        worker_started = deleteall_count > 0;
      });
  (void)worker_subscription;

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);
  ASSERT_TRUE(cacheMapAndWait(
      executor, node, makeLongRunningActionMap(),
      spawned.node_namespace))
    << "The shutdown test map was not admitted";
  auto goal_future = client->async_send_goal(makeLongRunningActionGoal());
  ASSERT_TRUE(waitForFuture(executor, goal_future, 10s))
    << "The shutdown test Action goal was not acknowledged";
  ASSERT_NE(goal_future.get(), nullptr)
    << "The shutdown test Action goal was rejected";

  // The Action handle proves admission; the first replacement marker proves
  // the owned worker has entered executePlanning().  Shutdown is therefore
  // synchronized to an active worker rather than an arbitrary delay.
  const auto worker_deadline = std::chrono::steady_clock::now() + 5s;
  while (!worker_started && std::chrono::steady_clock::now() < worker_deadline)
  {
    executor.spin_some(20ms);
    std::this_thread::sleep_for(1ms);
  }
  ASSERT_TRUE(worker_started)
    << "The Action worker did not enter planning before shutdown";

  ASSERT_EQ(kill(shutdown_pid, SIGINT), 0);

  int child_status = 0;
  pid_t wait_result = 0;
  bool wait_failed = false;
  const auto exit_deadline = std::chrono::steady_clock::now() + 10s;
  while (std::chrono::steady_clock::now() < exit_deadline)
  {
    wait_result = waitpid(shutdown_pid, &child_status, WNOHANG);
    if (wait_result == shutdown_pid)
      break;
    if (wait_result < 0)
    {
      if (errno == EINTR)
        continue;
      wait_failed = true;
      break;
    }
    executor.spin_some(20ms);
    std::this_thread::sleep_for(5ms);
  }

  if (wait_result == 0)
    wait_result = waitpid(shutdown_pid, &child_status, WNOHANG);

  if (wait_result == shutdown_pid)
    shutdown_process.disarm();
  EXPECT_FALSE(wait_failed) << "waitpid failed while waiting for graceful exit";
  ASSERT_EQ(wait_result, shutdown_pid)
    << "raystar_node did not exit within 10 seconds after SIGINT";
  ASSERT_TRUE(WIFEXITED(child_status));
  EXPECT_EQ(WEXITSTATUS(child_status), 0);

  const auto clear_deadline = std::chrono::steady_clock::now() + 1s;
  while (deleteall_count < 2 &&
      std::chrono::steady_clock::now() < clear_deadline)
  {
    executor.spin_some(20ms);
    std::this_thread::sleep_for(1ms);
  }
  // The first DELETEALL is the replacement marker that proves the worker was
  // active.  A second snapshot is published by the child destructor, but DDS
  // may drop an in-flight transient-local sample as the process exits.  Keep
  // this cross-process assertion at the observable best-effort boundary.
  EXPECT_GE(deleteall_count, 1u)
    << "The active worker did not publish a DELETEALL snapshot";
}

TEST_F(IntegrationTestFixture, PreservesContinuousEndpointsAndMetricCost)
{
  auto node = std::make_shared<rclcpp::Node>("test_continuous_endpoint_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto request = makeTestRequest();
  request->map.info.resolution = 0.25f;
  request->map.info.origin.position.x = -4.0;
  request->map.info.origin.position.y = 7.0;
  const auto to_world = [&](double gx, double gy) {
    return std::pair<double, double>{
      request->map.info.origin.position.x + gx * request->map.info.resolution,
      request->map.info.origin.position.y + gy * request->map.info.resolution};
  };
  const auto start_world = to_world(5.25, 15.25);
  const auto goal_world = to_world(25.75, 15.75);
  request->start.pose.position.x = start_world.first;
  request->start.pose.position.y = start_world.second;
  request->goal.pose.position.x = goal_world.first;
  request->goal.pose.position.y = goal_world.second;
  request->k = 1;

  const auto response = callService(executor, client, request);
  ASSERT_NE(response, nullptr);
  ASSERT_TRUE(response->success) << response->message;
  ASSERT_FALSE(response->path_results.empty());
  const auto& path_result = response->path_results.front();
  const auto& path = path_result.path;
  ASSERT_GE(path.poses.size(), 2u);
  EXPECT_NEAR(path.poses.front().pose.position.x, start_world.first, 1e-9);
  EXPECT_NEAR(path.poses.front().pose.position.y, start_world.second, 1e-9);
  EXPECT_NEAR(path.poses.back().pose.position.x, goal_world.first, 1e-9);
  EXPECT_NEAR(path.poses.back().pose.position.y, goal_world.second, 1e-9);

  double sampled_length = 0.0;
  for (size_t i = 1; i < path.poses.size(); ++i)
  {
    const auto& previous = path.poses[i - 1].pose.position;
    const auto& current = path.poses[i].pose.position;
    sampled_length += std::hypot(
      current.x - previous.x, current.y - previous.y);
  }
  EXPECT_NEAR(path_result.cost, sampled_length, 1e-8);
  EXPECT_GE(path_result.cost,
    std::hypot(goal_world.first - start_world.first,
               goal_world.second - start_world.second));
  expectIndependentCollisionFreePath(
    request->map, path, request->allow_unknown);
  EXPECT_EQ(response->result_info.status, PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_TRUE(response->result_info.request_satisfied);
  expectStructuredResultInvariants(*response);
}

TEST_F(IntegrationTestFixture,
  RejectsOccupiedAndForcedBorderEndpointsSymmetricallyThenRecovers)
{
  auto node = std::make_shared<rclcpp::Node>(
    "test_endpoint_policy_matrix_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  using Mutator = std::function<void(RaystarService::Request&)>;
  struct InvalidEndpointCase
  {
    const char* name;
    const char* expected_endpoint;
    const char* expected_reason;
    Mutator mutate;
  };
  const std::vector<InvalidEndpointCase> invalid_cases = {
    {"occupied start", "invalid start", "occupied", [](auto& request) {
      request.start.pose.position.x = 15.5;
      request.start.pose.position.y = 15.5;
    }},
    {"occupied goal", "invalid goal", "occupied", [](auto& request) {
      request.goal.pose.position.x = 15.5;
      request.goal.pose.position.y = 15.5;
    }},
    {"forced-border start", "invalid start", "map boundary", [](auto& request) {
      request.start.pose.position.x = 0.5;
      request.start.pose.position.y = 15.5;
    }},
    {"forced-border goal", "invalid goal", "map boundary", [](auto& request) {
      request.goal.pose.position.x = 29.5;
      request.goal.pose.position.y = 15.5;
    }},
  };

  for (const auto& invalid_case : invalid_cases)
  {
    SCOPED_TRACE(invalid_case.name);
    auto request = makeTestRequest();
    request->k = 1;
    invalid_case.mutate(*request);
    const auto response = callService(executor, client, request);

    ASSERT_NE(response, nullptr);
    EXPECT_FALSE(response->success) << response->message;
    EXPECT_TRUE(response->path_results.empty());
    EXPECT_TRUE(response->debug_nodes.empty());
    EXPECT_EQ(response->result_info.status,
      PlanningResultInfo::STATUS_INVALID_REQUEST);
    EXPECT_EQ(response->result_info.limits_reached,
      PlanningResultInfo::LIMIT_NONE);
    EXPECT_EQ(response->result_info.requested_path_count, 1u);
    EXPECT_EQ(response->result_info.found_path_count, 0u);
    EXPECT_EQ(response->result_info.returned_path_count, 0u);
    EXPECT_EQ(response->result_info.expanded_nodes, 0u);
    EXPECT_FALSE(response->result_info.request_satisfied);
    EXPECT_FALSE(response->result_info.search_complete);
    expectStructuredResultInvariants(*response);
    const auto diagnostic = lowercase(response->message);
    EXPECT_NE(diagnostic.find(invalid_case.expected_endpoint),
      std::string::npos) << response->message;
    EXPECT_NE(diagnostic.find(invalid_case.expected_reason),
      std::string::npos) << response->message;
  }

  auto recovery_request = makeTestRequest();
  recovery_request->k = 1;
  const auto recovered = callService(executor, client, recovery_request);
  ASSERT_NE(recovered, nullptr);
  ASSERT_TRUE(recovered->success) << recovered->message;
  ASSERT_EQ(recovered->path_results.size(), 1u);
  EXPECT_EQ(recovered->result_info.status,
    PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_EQ(recovered->result_info.limits_reached,
    PlanningResultInfo::LIMIT_NONE);
  EXPECT_EQ(recovered->result_info.requested_path_count, 1u);
  EXPECT_EQ(recovered->result_info.found_path_count, 1u);
  EXPECT_EQ(recovered->result_info.returned_path_count, 1u);
  EXPECT_TRUE(recovered->result_info.request_satisfied);
  EXPECT_TRUE(recovered->result_info.search_complete);
  EXPECT_TRUE(recovered->result_info.output_complete);
  expectStructuredResultInvariants(*recovered);
}

TEST_F(IntegrationTestFixture, RejectedAndNoPathRequestsDoNotPoisonServer)
{
  auto node = std::make_shared<rclcpp::Node>("test_invalid_map_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));

  ASSERT_TRUE(client->wait_for_service(2s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto call_service = [&](const RaystarService::Request::SharedPtr& request)
    -> RaystarService::Response::SharedPtr
  {
    auto future = client->async_send_request(request);
    const auto end_time = std::chrono::steady_clock::now() + 10s;
    while (std::chrono::steady_clock::now() < end_time)
    {
      if (future.wait_for(100ms) == std::future_status::ready) break;
      executor.spin_some(100ms);
    }
    if (future.wait_for(0s) != std::future_status::ready)
    {
      ADD_FAILURE() << "Service call did not complete within 10 seconds";
      return nullptr;
    }
    return future.get();
  };

  auto short_request = makeTestRequest();
  short_request->map.data.pop_back();
  auto short_result = call_service(short_request);
  ASSERT_NE(short_result, nullptr);
  EXPECT_FALSE(short_result->success);
  EXPECT_NE(short_result->message.find("data size"), std::string::npos);

  auto long_request = makeTestRequest();
  long_request->map.data.push_back(0);
  auto long_result = call_service(long_request);
  ASSERT_NE(long_result, nullptr);
  EXPECT_FALSE(long_result->success);
  EXPECT_NE(long_result->message.find("data size"), std::string::npos);

  auto invalid_negative_request = makeTestRequest();
  invalid_negative_request->allow_unknown = true;
  invalid_negative_request->map.data[0] = static_cast<int8_t>(-2);
  auto invalid_negative_result = call_service(invalid_negative_request);
  ASSERT_NE(invalid_negative_result, nullptr);
  EXPECT_FALSE(invalid_negative_result->success);
  EXPECT_NE(invalid_negative_result->message.find("occupancy value -2"),
    std::string::npos) << invalid_negative_result->message;
  EXPECT_NE(invalid_negative_result->message.find("[-1, 100]"),
    std::string::npos) << invalid_negative_result->message;

  auto invalid_high_request = makeTestRequest();
  invalid_high_request->map.data[0] = static_cast<int8_t>(101);
  auto invalid_high_result = call_service(invalid_high_request);
  ASSERT_NE(invalid_high_result, nullptr);
  EXPECT_FALSE(invalid_high_result->success);
  EXPECT_NE(invalid_high_result->message.find("occupancy value 101"),
    std::string::npos) << invalid_high_result->message;
  EXPECT_NE(invalid_high_result->message.find("[-1, 100]"),
    std::string::npos) << invalid_high_result->message;

  auto blocked_request = makeTestRequest();
  for (unsigned int y = 0; y < blocked_request->map.info.height; ++y)
    blocked_request->map.data[y * blocked_request->map.info.width + 15] = 100;
  auto blocked_result = call_service(blocked_request);
  ASSERT_NE(blocked_result, nullptr);
  EXPECT_FALSE(blocked_result->success);
  EXPECT_NE(blocked_result->message.find("No path exists"), std::string::npos);
  EXPECT_TRUE(blocked_result->debug_nodes.empty());
  EXPECT_EQ(blocked_result->result_info.status,
    PlanningResultInfo::STATUS_NO_PATH);
  EXPECT_TRUE(blocked_result->result_info.search_complete);
  EXPECT_TRUE(blocked_result->result_info.output_complete);
  EXPECT_FALSE(blocked_result->result_info.request_satisfied);
  EXPECT_EQ(blocked_result->result_info.requested_path_count, 3u);
  EXPECT_EQ(blocked_result->result_info.found_path_count, 0u);
  expectStructuredResultInvariants(*blocked_result);

  auto shared_vertex_request = makeTestRequest();
  std::fill(shared_vertex_request->map.data.begin(),
    shared_vertex_request->map.data.end(), 0);
  shared_vertex_request->map.data[10 * shared_vertex_request->map.info.width + 10] = 100;
  shared_vertex_request->map.data[11 * shared_vertex_request->map.info.width + 11] = 100;
  auto shared_vertex_result = call_service(shared_vertex_request);
  ASSERT_NE(shared_vertex_result, nullptr);
  EXPECT_FALSE(shared_vertex_result->success);
  EXPECT_NE(shared_vertex_result->message.find(
    "Unsupported obstacle topology"), std::string::npos)
    << shared_vertex_result->message;
  EXPECT_NE(shared_vertex_result->message.find("(11, 11)"), std::string::npos)
    << shared_vertex_result->message;
  EXPECT_TRUE(shared_vertex_result->debug_nodes.empty());
  EXPECT_EQ(shared_vertex_result->result_info.status,
    PlanningResultInfo::STATUS_FAILED);

  auto valid_result = call_service(makeTestRequest());
  ASSERT_NE(valid_result, nullptr);
  EXPECT_TRUE(valid_result->success);
  EXPECT_GE(valid_result->path_results.size(), 1u);
  expectStructuredResultInvariants(*valid_result);
}

TEST_F(IntegrationTestFixture, RejectsInvalidFramesMetricMetadataAndPosesThenRecovers)
{
  auto node = std::make_shared<rclcpp::Node>("test_geometry_contract_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  using Mutator = std::function<void(RaystarService::Request&)>;
  struct InvalidCase
  {
    const char* name;
    const char* expected_keyword;
    Mutator mutate;
  };

  const double nan = std::numeric_limits<double>::quiet_NaN();
  const double infinity = std::numeric_limits<double>::infinity();
  const float float_nan = std::numeric_limits<float>::quiet_NaN();
  const float float_infinity = std::numeric_limits<float>::infinity();
  constexpr double sqrt_half = 0.70710678118654752440;

  const std::vector<InvalidCase> invalid_cases = {
    {"empty map frame", "frame", [](auto& request) {
      request.map.header.frame_id.clear();
    }},
    {"empty start frame", "frame", [](auto& request) {
      request.start.header.frame_id.clear();
    }},
    {"empty goal frame", "frame", [](auto& request) {
      request.goal.header.frame_id.clear();
    }},
    {"mismatched start frame", "frame", [](auto& request) {
      request.start.header.frame_id = "odom";
    }},
    {"mismatched goal frame", "frame", [](auto& request) {
      request.goal.header.frame_id = "odom";
    }},
    {"oversized map frame", "256", [](auto& request) {
      request.map.header.frame_id.assign(257, 'm');
    }},
    {"zero resolution", "resolution", [](auto& request) {
      request.map.info.resolution = 0.0f;
    }},
    {"negative resolution", "resolution", [](auto& request) {
      request.map.info.resolution = -1.0f;
    }},
    {"NaN resolution", "resolution", [float_nan](auto& request) {
      request.map.info.resolution = float_nan;
    }},
    {"infinite resolution", "resolution", [float_infinity](auto& request) {
      request.map.info.resolution = float_infinity;
    }},
    {"NaN origin x", "origin", [nan](auto& request) {
      request.map.info.origin.position.x = nan;
    }},
    {"infinite origin y", "origin", [infinity](auto& request) {
      request.map.info.origin.position.y = infinity;
    }},
    {"nonzero origin z", "origin", [](auto& request) {
      request.map.info.origin.position.z = 1.0;
    }},
    {"NaN origin z", "origin", [nan](auto& request) {
      request.map.info.origin.position.z = nan;
    }},
    {"zero origin quaternion", "orientation", [](auto& request) {
      request.map.info.origin.orientation.x = 0.0;
      request.map.info.origin.orientation.y = 0.0;
      request.map.info.origin.orientation.z = 0.0;
      request.map.info.origin.orientation.w = 0.0;
    }},
    {"non-unit origin quaternion", "orientation", [](auto& request) {
      request.map.info.origin.orientation.w = 2.0;
    }},
    {"nonzero map yaw", "orientation", [sqrt_half](auto& request) {
      request.map.info.origin.orientation.z = sqrt_half;
      request.map.info.origin.orientation.w = sqrt_half;
    }},
    {"nonzero map roll", "orientation", [sqrt_half](auto& request) {
      request.map.info.origin.orientation.x = sqrt_half;
      request.map.info.origin.orientation.w = sqrt_half;
    }},
    {"NaN map quaternion", "orientation", [nan](auto& request) {
      request.map.info.origin.orientation.z = nan;
    }},
    {"NaN start x", "start", [nan](auto& request) {
      request.start.pose.position.x = nan;
    }},
    {"infinite start y", "start", [infinity](auto& request) {
      request.start.pose.position.y = infinity;
    }},
    {"nonzero start z", "start", [](auto& request) {
      request.start.pose.position.z = 1.0;
    }},
    {"infinite goal x", "goal", [infinity](auto& request) {
      request.goal.pose.position.x = -infinity;
    }},
    {"NaN goal y", "goal", [nan](auto& request) {
      request.goal.pose.position.y = nan;
    }},
    {"nonzero goal z", "goal", [](auto& request) {
      request.goal.pose.position.z = -1.0;
    }},
  };

  for (const auto& invalid_case : invalid_cases)
  {
    SCOPED_TRACE(invalid_case.name);
    auto request = makeTestRequest();
    invalid_case.mutate(*request);

    auto response = callService(executor, client, request);
    ASSERT_NE(response, nullptr);
    EXPECT_FALSE(response->success) << response->message;
    EXPECT_TRUE(response->path_results.empty());
    EXPECT_TRUE(response->debug_nodes.empty());
    EXPECT_EQ(response->result_info.status,
      PlanningResultInfo::STATUS_INVALID_REQUEST);
    expectStructuredResultInvariants(*response);
    EXPECT_NE(lowercase(response->message).find(invalid_case.expected_keyword),
      std::string::npos) << response->message;
  }

  auto valid_response = callService(executor, client, makeTestRequest());
  ASSERT_NE(valid_response, nullptr);
  EXPECT_TRUE(valid_response->success) << valid_response->message;
  EXPECT_FALSE(valid_response->path_results.empty());
}

TEST_F(IntegrationTestFixture, RejectsKAboveConfiguredMaximumThenRecovers)
{
  auto node = std::make_shared<rclcpp::Node>("test_max_k_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto excessive_request = makeTestRequest();
  excessive_request->k = 101;
  auto excessive_response = callService(executor, client, excessive_request);
  ASSERT_NE(excessive_response, nullptr);
  EXPECT_FALSE(excessive_response->success);
  EXPECT_TRUE(excessive_response->path_results.empty());
  EXPECT_TRUE(excessive_response->debug_nodes.empty());
  EXPECT_EQ(excessive_response->result_info.status,
    PlanningResultInfo::STATUS_INVALID_REQUEST);
  EXPECT_NE(excessive_response->message.find("max_k=100"), std::string::npos)
    << excessive_response->message;

  auto valid_response = callService(executor, client, makeTestRequest());
  ASSERT_NE(valid_response, nullptr);
  EXPECT_TRUE(valid_response->success) << valid_response->message;
  EXPECT_FALSE(valid_response->path_results.empty());
}

TEST_F(IntegrationTestFixture, DynamicParametersAreDescribedAndUpdatedAtomically)
{
  const auto spawned = spawnRaystarNode("dynamic_parameters");
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec dynamic-parameters raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard parameter_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_dynamic_parameter_client");
  auto service_client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  auto parameter_client = std::make_shared<rclcpp::AsyncParametersClient>(
    node, raystarEndpoint(spawned.node_namespace));
  ASSERT_TRUE(service_client->wait_for_service(5s));
  ASSERT_TRUE(parameter_client->wait_for_service(5s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  struct ExpectedDescriptor
  {
    const char* name;
    int64_t minimum;
    int64_t maximum;
    bool read_only;
  };
  const int64_t max_size_t_parameter = []() constexpr {
    if constexpr (sizeof(size_t) < sizeof(int64_t))
      return static_cast<int64_t>(std::numeric_limits<size_t>::max());
    return std::numeric_limits<int64_t>::max();
  }();
  const int64_t max_int =
    static_cast<int64_t>(std::numeric_limits<int>::max());
  const int64_t max_timeout = std::chrono::milliseconds::max().count() - 1;
  const int64_t max_timer_period =
    std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::nanoseconds::max() - std::chrono::milliseconds(1)).count();
  const std::vector<ExpectedDescriptor> expected_descriptors = {
    {"occupied_threshold", 1, 100, false},
    {"max_k", 1, max_int, false},
    {"max_nodes", 1, max_int, false},
    {"planning_timeout_ms", 1, max_timeout, false},
    {"max_map_cells", 1, max_int, false},
    {"max_map_bytes", 32, max_size_t_parameter, false},
    {"max_path_points", 2, max_int, false},
    {"max_debug_nodes", 0, max_int, false},
    {"max_response_bytes", 1024, max_size_t_parameter, false},
    {"path_visualization_republish_period_ms", 0, max_timer_period, true},
  };
  std::vector<std::string> parameter_names;
  parameter_names.reserve(expected_descriptors.size());
  for (const auto& expected : expected_descriptors)
    parameter_names.emplace_back(expected.name);

  auto descriptor_future = parameter_client->describe_parameters(parameter_names);
  ASSERT_TRUE(waitForFuture(executor, descriptor_future, 5s));
  const auto descriptors = descriptor_future.get();
  ASSERT_EQ(descriptors.size(), expected_descriptors.size());
  for (size_t i = 0; i < descriptors.size(); ++i)
  {
    SCOPED_TRACE(expected_descriptors[i].name);
    const auto& descriptor = descriptors[i];
    EXPECT_EQ(descriptor.name, expected_descriptors[i].name);
    EXPECT_EQ(descriptor.type,
      static_cast<uint8_t>(rclcpp::ParameterType::PARAMETER_INTEGER));
    EXPECT_FALSE(descriptor.description.empty());
    EXPECT_EQ(descriptor.read_only, expected_descriptors[i].read_only);
    ASSERT_EQ(descriptor.integer_range.size(), 1u);
    EXPECT_EQ(descriptor.integer_range.front().from_value,
      expected_descriptors[i].minimum);
    EXPECT_EQ(descriptor.integer_range.front().to_value,
      expected_descriptors[i].maximum);
    EXPECT_EQ(descriptor.integer_range.front().step, 1u);
  }

  auto initial_future = parameter_client->get_parameters(
    {"max_k", "occupied_threshold",
     "path_visualization_republish_period_ms"});
  ASSERT_TRUE(waitForFuture(executor, initial_future, 5s));
  const auto initial_values = initial_future.get();
  ASSERT_EQ(initial_values.size(), 3u);
  EXPECT_EQ(initial_values[0].as_int(), 100);
  EXPECT_EQ(initial_values[1].as_int(), 99);
  EXPECT_EQ(initial_values[2].as_int(), 2000);

  auto valid_update_future = parameter_client->set_parameters_atomically({
    rclcpp::Parameter("max_k", 1),
    rclcpp::Parameter("occupied_threshold", 50),
  });
  ASSERT_TRUE(waitForFuture(executor, valid_update_future, 5s));
  const auto valid_update = valid_update_future.get();
  ASSERT_TRUE(valid_update.successful) << valid_update.reason;

  auto updated_values_future = parameter_client->get_parameters(
    {"max_k", "occupied_threshold"});
  ASSERT_TRUE(waitForFuture(executor, updated_values_future, 5s));
  const auto updated_values = updated_values_future.get();
  ASSERT_EQ(updated_values.size(), 2u);
  EXPECT_EQ(updated_values[0].as_int(), 1);
  EXPECT_EQ(updated_values[1].as_int(), 50);

  auto excessive_k_request = makeTestRequest();
  excessive_k_request->k = 2;
  auto excessive_k_response = callService(
    executor, service_client, excessive_k_request);
  ASSERT_NE(excessive_k_response, nullptr);
  EXPECT_FALSE(excessive_k_response->success);
  EXPECT_NE(excessive_k_response->message.find("max_k=1"), std::string::npos)
    << excessive_k_response->message;

  auto below_threshold_response = callService(
    executor, service_client, makeVerticalBarrierRequest(49));
  ASSERT_NE(below_threshold_response, nullptr);
  EXPECT_TRUE(below_threshold_response->success)
    << below_threshold_response->message;
  auto at_threshold_response = callService(
    executor, service_client, makeVerticalBarrierRequest(50));
  ASSERT_NE(at_threshold_response, nullptr);
  EXPECT_FALSE(at_threshold_response->success);
  EXPECT_EQ(at_threshold_response->result_info.status,
    PlanningResultInfo::STATUS_NO_PATH);

  auto invalid_update_future = parameter_client->set_parameters_atomically({
    rclcpp::Parameter("max_k", 0),
  });
  ASSERT_TRUE(waitForFuture(executor, invalid_update_future, 5s));
  const auto invalid_update = invalid_update_future.get();
  EXPECT_FALSE(invalid_update.successful);
  EXPECT_FALSE(invalid_update.reason.empty());
  auto after_invalid_future = parameter_client->get_parameters({"max_k"});
  ASSERT_TRUE(waitForFuture(executor, after_invalid_future, 5s));
  const auto after_invalid = after_invalid_future.get();
  ASSERT_EQ(after_invalid.size(), 1u);
  EXPECT_EQ(after_invalid.front().as_int(), 1);

  auto mixed_update_future = parameter_client->set_parameters_atomically({
    rclcpp::Parameter("max_k", 2),
    rclcpp::Parameter("occupied_threshold", 101),
  });
  ASSERT_TRUE(waitForFuture(executor, mixed_update_future, 5s));
  const auto mixed_update = mixed_update_future.get();
  EXPECT_FALSE(mixed_update.successful);
  EXPECT_FALSE(mixed_update.reason.empty());

  auto unchanged_values_future = parameter_client->get_parameters(
    {"max_k", "occupied_threshold"});
  ASSERT_TRUE(waitForFuture(executor, unchanged_values_future, 5s));
  const auto unchanged_values = unchanged_values_future.get();
  ASSERT_EQ(unchanged_values.size(), 2u);
  EXPECT_EQ(unchanged_values[0].as_int(), 1);
  EXPECT_EQ(unchanged_values[1].as_int(), 50);

  auto read_only_update_future = parameter_client->set_parameters_atomically({
    rclcpp::Parameter("path_visualization_republish_period_ms", 0),
  });
  ASSERT_TRUE(waitForFuture(executor, read_only_update_future, 5s));
  const auto read_only_update = read_only_update_future.get();
  EXPECT_FALSE(read_only_update.successful);
  EXPECT_FALSE(read_only_update.reason.empty());
  auto period_future = parameter_client->get_parameters(
    {"path_visualization_republish_period_ms"});
  ASSERT_TRUE(waitForFuture(executor, period_future, 5s));
  const auto period = period_future.get();
  ASSERT_EQ(period.size(), 1u);
  EXPECT_EQ(period.front().as_int(), 2000);

  auto recovery_update_future = parameter_client->set_parameters_atomically({
    rclcpp::Parameter("max_k", 100),
    rclcpp::Parameter("occupied_threshold", 99),
  });
  ASSERT_TRUE(waitForFuture(executor, recovery_update_future, 5s));
  const auto recovery_update = recovery_update_future.get();
  ASSERT_TRUE(recovery_update.successful) << recovery_update.reason;
  auto recovered_values_future = parameter_client->get_parameters(
    {"max_k", "occupied_threshold"});
  ASSERT_TRUE(waitForFuture(executor, recovered_values_future, 5s));
  const auto recovered_values = recovered_values_future.get();
  ASSERT_EQ(recovered_values.size(), 2u);
  EXPECT_EQ(recovered_values[0].as_int(), 100);
  EXPECT_EQ(recovered_values[1].as_int(), 99);
  auto recovery_response = callService(
    executor, service_client, makeTestRequest());
  ASSERT_NE(recovery_response, nullptr);
  EXPECT_TRUE(recovery_response->success) << recovery_response->message;
  expectStructuredResultInvariants(*recovery_response);
}

TEST_F(IntegrationTestFixture, ConfiguredOccupancyThresholdIsApplied)
{
  const auto spawned = spawnRaystarNode(
    "occupancy_threshold_50", {"occupied_threshold:=50"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec occupancy-threshold raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard threshold_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>(
    "test_configured_occupancy_threshold_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto below_threshold = callService(
    executor, client, makeVerticalBarrierRequest(49));
  ASSERT_NE(below_threshold, nullptr);
  EXPECT_TRUE(below_threshold->success) << below_threshold->message;

  auto at_threshold = callService(
    executor, client, makeVerticalBarrierRequest(50));
  ASSERT_NE(at_threshold, nullptr);
  EXPECT_FALSE(at_threshold->success);
  EXPECT_NE(at_threshold->message.find("No path exists"), std::string::npos)
    << at_threshold->message;
}

TEST_F(IntegrationTestFixture, RejectsInvalidOccupancyThresholdAtStartup)
{
  struct InvalidThresholdCase
  {
    const char* namespace_remap;
    const char* parameter_override;
  };
  const std::vector<InvalidThresholdCase> invalid_cases = {
    {"occupancy_threshold_zero", "occupied_threshold:=0"},
    {"occupancy_threshold_101", "occupied_threshold:=101"},
  };

  for (const auto& invalid_case : invalid_cases)
  {
    SCOPED_TRACE(invalid_case.parameter_override);
    const auto spawned = spawnRaystarNode(
      invalid_case.namespace_remap, {invalid_case.parameter_override});
    ASSERT_GT(spawned.pid, 0)
      << "Unable to exec invalid-parameter raystar_node: "
      << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
      "fork/pipe failure");
    ChildProcessGuard invalid_process(spawned.pid);
    int status = 0;
    ASSERT_TRUE(waitForChildExit(spawned.pid, status, 5s))
      << "raystar_node accepted " << invalid_case.parameter_override;
    invalid_process.disarm();
    ASSERT_TRUE(WIFEXITED(status));
    EXPECT_EQ(WEXITSTATUS(status), 1);
  }
}

TEST_F(IntegrationTestFixture, MaxNodesStopsExpansionAndLimitedServerRecovers)
{
  const auto spawned = spawnRaystarNode(
    "resource_limited", {"max_nodes:=1", "planning_timeout_ms:=60000",
      "max_debug_nodes:=1"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec max-nodes raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard limited_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_max_nodes_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto limited_request = makeTestRequest();
  limited_request->include_debug = true;
  auto limited_response = callService(executor, client, limited_request);
  ASSERT_NE(limited_response, nullptr);
  EXPECT_FALSE(limited_response->success);
  EXPECT_TRUE(limited_response->path_results.empty());
  EXPECT_EQ(limited_response->debug_nodes.size(), 1u);
  EXPECT_EQ(limited_response->result_info.status,
    PlanningResultInfo::STATUS_PARTIAL_SEARCH);
  EXPECT_NE(limited_response->result_info.limits_reached &
    PlanningResultInfo::LIMIT_MAX_NODES, 0u);
  EXPECT_FALSE(limited_response->result_info.search_complete);
  EXPECT_TRUE(limited_response->result_info.output_complete);
  EXPECT_FALSE(limited_response->result_info.request_satisfied);
  EXPECT_EQ(limited_response->result_info.requested_path_count, 3u);
  EXPECT_EQ(limited_response->result_info.found_path_count, 0u);
  EXPECT_EQ(limited_response->result_info.returned_path_count, 0u);
  expectStructuredResultInvariants(*limited_response);
  EXPECT_NE(limited_response->message.find("max_nodes=1"), std::string::npos)
    << limited_response->message;

  auto root_solution_request = makeTestRequest();
  std::fill(root_solution_request->map.data.begin(),
    root_solution_request->map.data.end(), 0);
  root_solution_request->k = 1;
  root_solution_request->include_debug = true;
  auto recovered_response = callService(
    executor, client, root_solution_request);
  ASSERT_NE(recovered_response, nullptr);
  EXPECT_TRUE(recovered_response->success) << recovered_response->message;
  EXPECT_EQ(recovered_response->path_results.size(), 1u);
  EXPECT_EQ(recovered_response->debug_nodes.size(), 1u);
  EXPECT_EQ(recovered_response->result_info.status,
    PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_TRUE(recovered_response->result_info.request_satisfied);
  expectStructuredResultInvariants(*recovered_response);
}

TEST_F(IntegrationTestFixture, PlanningTimeoutHasStructuredLimitStatus)
{
  const auto spawned = spawnRaystarNode(
    "planning_timeout", {"planning_timeout_ms:=1"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec planning-timeout raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard timeout_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_planning_timeout_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto goal = makeLongRunningActionGoal();
  auto request = std::make_shared<RaystarService::Request>();
  request->map = makeLongRunningActionMap();
  request->start = std::move(goal.start);
  request->goal = std::move(goal.goal);
  request->k = goal.k;
  request->allow_self_crossing = goal.allow_self_crossing;
  request->allow_unknown = goal.allow_unknown;
  request->include_debug = false;

  const auto response = callService(executor, client, request);
  ASSERT_NE(response, nullptr);
  EXPECT_FALSE(response->success);
  EXPECT_EQ(response->result_info.status,
    PlanningResultInfo::STATUS_PARTIAL_SEARCH);
  EXPECT_NE(response->result_info.limits_reached &
    PlanningResultInfo::LIMIT_TIMEOUT, 0u);
  EXPECT_FALSE(response->result_info.search_complete);
  EXPECT_FALSE(response->result_info.request_satisfied);
  EXPECT_TRUE(response->path_results.empty());
  expectStructuredResultInvariants(*response);
}

TEST_F(IntegrationTestFixture, MapCellAdmissionHappensBeforePlannerAndRecovers)
{
  const auto spawned = spawnRaystarNode(
    "map_cell_limited", {"max_map_cells:=899"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec map-cell-limit raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard limited_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_map_cell_limit_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto excessive = makeTestRequest();  // 30 * 30 = 900 cells
  auto rejected = callService(executor, client, excessive);
  ASSERT_NE(rejected, nullptr);
  EXPECT_FALSE(rejected->success);
  EXPECT_TRUE(rejected->path_results.empty());
  EXPECT_NE(rejected->message.find("max_map_cells=899"), std::string::npos)
    << rejected->message;

  auto smaller = makeTestRequest();
  smaller->map.info.width = 29;
  smaller->map.data.assign(29 * 30, 0);
  for (unsigned int y = 10; y < 20; ++y)
    for (unsigned int x = 10; x < 20; ++x)
      smaller->map.data[y * 29 + x] = 100;
  auto recovered = callService(executor, client, smaller);
  ASSERT_NE(recovered, nullptr);
  EXPECT_TRUE(recovered->success) << recovered->message;
  EXPECT_FALSE(recovered->path_results.empty());
}

TEST_F(IntegrationTestFixture, MapByteAdmissionUsesWorkingSetEstimate)
{
  const auto spawned = spawnRaystarNode(
    "map_byte_limited", {"max_map_cells:=1000", "max_map_bytes:=28000"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec map-byte-limit raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard limited_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_map_byte_limit_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto excessive = makeTestRequest();
  auto rejected = callService(executor, client, excessive);
  ASSERT_NE(rejected, nullptr);
  EXPECT_FALSE(rejected->success);
  EXPECT_NE(rejected->message.find("max_map_bytes=28000"), std::string::npos)
    << rejected->message;

  auto smaller = makeTestRequest();  // 29 * 30 * 32 = 27840 bytes
  smaller->map.info.width = 29;
  smaller->map.data.assign(29 * 30, 0);
  for (unsigned int y = 10; y < 20; ++y)
    for (unsigned int x = 10; x < 20; ++x)
      smaller->map.data[y * 29 + x] = 100;
  auto recovered = callService(executor, client, smaller);
  ASSERT_NE(recovered, nullptr);
  EXPECT_TRUE(recovered->success) << recovered->message;
}

TEST_F(IntegrationTestFixture, PathInterpolationLimitRejectsOversizedPath)
{
  const auto spawned = spawnRaystarNode(
    "path_points_limited", {"max_path_points:=10"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec path-point-limit raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard limited_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_path_points_limit_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto excessive = makeTestRequest();
  excessive->k = 1;
  auto rejected = callService(executor, client, excessive);
  ASSERT_NE(rejected, nullptr);
  EXPECT_FALSE(rejected->success);
  EXPECT_TRUE(rejected->path_results.empty());
  EXPECT_EQ(rejected->result_info.status,
    PlanningResultInfo::STATUS_PARTIAL_OUTPUT);
  EXPECT_NE(rejected->result_info.limits_reached &
    PlanningResultInfo::LIMIT_MAX_PATH_POINTS, 0u);
  EXPECT_TRUE(rejected->result_info.search_complete);
  EXPECT_FALSE(rejected->result_info.output_complete);
  EXPECT_EQ(rejected->result_info.requested_path_count, 1u);
  EXPECT_EQ(rejected->result_info.found_path_count, 1u);
  EXPECT_EQ(rejected->result_info.returned_path_count, 0u);
  expectStructuredResultInvariants(*rejected);
  EXPECT_NE(rejected->message.find("max_path_points=10"), std::string::npos)
    << rejected->message;

  auto short_request = makeTestRequest();
  std::fill(short_request->map.data.begin(), short_request->map.data.end(), 0);
  short_request->k = 1;
  short_request->goal.pose.position.x = 6.5;
  auto recovered = callService(executor, client, short_request);
  ASSERT_NE(recovered, nullptr);
  EXPECT_TRUE(recovered->success) << recovered->message;
  ASSERT_EQ(recovered->path_results.size(), 1u);
  EXPECT_LE(recovered->path_results.front().path.poses.size(), 10u);
}

TEST_F(IntegrationTestFixture, PathPointLimitIsCumulativeAcrossResponse)
{
  const auto spawned = spawnRaystarNode(
    "total_path_points_limited", {"max_path_points:=30"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec cumulative path-point-limit raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard limited_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>(
    "test_total_path_points_limit_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto request = makeTestRequest();
  request->k = 2;
  auto response = callService(executor, client, request);
  ASSERT_NE(response, nullptr);
  ASSERT_TRUE(response->success) << response->message;
  ASSERT_EQ(response->path_results.size(), 1u) << response->message;
  size_t total_points = 0;
  for (const auto& path_result : response->path_results)
    total_points += path_result.path.poses.size();
  EXPECT_LE(total_points, 30u);
  EXPECT_EQ(response->result_info.status,
    PlanningResultInfo::STATUS_PARTIAL_OUTPUT);
  EXPECT_NE(response->result_info.limits_reached &
    PlanningResultInfo::LIMIT_MAX_PATH_POINTS, 0u);
  EXPECT_TRUE(response->result_info.search_complete);
  EXPECT_FALSE(response->result_info.output_complete);
  EXPECT_FALSE(response->result_info.request_satisfied);
  EXPECT_EQ(response->result_info.requested_path_count, 2u);
  EXPECT_EQ(response->result_info.found_path_count, 2u);
  EXPECT_EQ(response->result_info.returned_path_count, 1u);
  expectStructuredResultInvariants(*response);
  EXPECT_NE(response->message.find("per-response max_path_points=30"),
    std::string::npos) << response->message;
}

TEST_F(IntegrationTestFixture, DebugNodeLimitBoundsResponseArrays)
{
  const auto spawned = spawnRaystarNode(
    "debug_nodes_limited", {"max_debug_nodes:=0"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec debug-node-limit raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard limited_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_debug_nodes_limit_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto request = makeTestRequest();
  request->include_debug = true;
  auto response = callService(executor, client, request);
  ASSERT_NE(response, nullptr);
  EXPECT_TRUE(response->success) << response->message;
  EXPECT_TRUE(response->debug_nodes.empty());
  EXPECT_TRUE(response->result_info.debug_requested);
  EXPECT_FALSE(response->result_info.debug_output_complete);
  EXPECT_GT(response->result_info.expanded_nodes, 0u);
  expectStructuredResultInvariants(*response);
  EXPECT_NE(response->message.find("Debug output limited"), std::string::npos)
    << response->message;
}

TEST_F(IntegrationTestFixture, ResponseByteLimitBoundsPathsAndRecovers)
{
  const auto spawned = spawnRaystarNode(
    "response_bytes_limited", {"max_response_bytes:=1024"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec response-byte-limit raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard limited_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_response_bytes_limit_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto excessive = makeTestRequest();
  excessive->k = 1;
  auto rejected = callService(executor, client, excessive);
  ASSERT_NE(rejected, nullptr);
  EXPECT_FALSE(rejected->success);
  EXPECT_TRUE(rejected->path_results.empty());
  EXPECT_EQ(rejected->result_info.status,
    PlanningResultInfo::STATUS_PARTIAL_OUTPUT);
  EXPECT_NE(rejected->result_info.limits_reached &
    PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES, 0u);
  EXPECT_TRUE(rejected->result_info.search_complete);
  EXPECT_FALSE(rejected->result_info.output_complete);
  EXPECT_EQ(rejected->result_info.requested_path_count, 1u);
  EXPECT_EQ(rejected->result_info.found_path_count, 1u);
  EXPECT_EQ(rejected->result_info.returned_path_count, 0u);
  expectStructuredResultInvariants(*rejected);
  EXPECT_NE(rejected->message.find("max_response_bytes=1024"), std::string::npos)
    << rejected->message;

  auto short_request = makeTestRequest();
  std::fill(short_request->map.data.begin(), short_request->map.data.end(), 0);
  short_request->k = 1;
  short_request->goal.pose.position.x = 6.5;
  auto recovered = callService(executor, client, short_request);
  ASSERT_NE(recovered, nullptr);
  EXPECT_TRUE(recovered->success) << recovered->message;
  EXPECT_EQ(recovered->path_results.size(), 1u);
}

TEST_F(IntegrationTestFixture, ReservedNoDeadlineTimeoutFailsAtStartup)
{
  const auto spawned = spawnRaystarNode(
    "invalid_timeout", {"planning_timeout_ms:=9223372036854775807"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec invalid-timeout raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard invalid_timeout_process(spawned.pid);
  int status = 0;
  ASSERT_TRUE(waitForChildExit(spawned.pid, status, 5s));
  invalid_timeout_process.disarm();
  ASSERT_TRUE(WIFEXITED(status));
  EXPECT_EQ(WEXITSTATUS(status), 1);
}

TEST_F(IntegrationTestFixture, MarkerSnapshotsReplaceClearAndRemainDurable)
{
  const auto spawned = spawnRaystarNode(
    "marker_lifecycle", {"planning_timeout_ms:=60000",
      "path_visualization_republish_period_ms:=200", "max_debug_nodes:=100"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec marker-lifecycle raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard marker_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_marker_lifecycle_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  const int64_t first_request_start_ns = node->now().nanoseconds();
  auto first_request = makeTestRequest();
  first_request->k = 3;
  // Geometry/CDT/tree marker topics are explicit opt-in debug output.
  first_request->include_debug = true;
  auto first_response = callService(executor, client, first_request);
  ASSERT_NE(first_response, nullptr);
  ASSERT_TRUE(first_response->success) << first_response->message;
  ASSERT_GT(first_response->path_results.size(), 1u);
  for (size_t path_index = 0;
    path_index < first_response->path_results.size(); ++path_index)
  {
    SCOPED_TRACE("returned path " + std::to_string(path_index));
    expectIndependentCollisionFreePath(
      first_request->map,
      first_response->path_results[path_index].path,
      first_request->allow_unknown);
  }

  using MarkerArray = visualization_msgs::msg::MarkerArray;
  using MarkerHistory = std::vector<MarkerArray::ConstSharedPtr>;
  MarkerHistory path_messages;
  MarkerHistory obstacle_messages;
  MarkerHistory cdt_messages;
  MarkerHistory tree_messages;

  const int64_t subscription_start_ns = node->now().nanoseconds();
  std::vector<rclcpp::SubscriptionBase::SharedPtr> subscriptions;
  auto marker_qos = rclcpp::QoS(rclcpp::KeepLast(10)).transient_local();
  subscriptions.push_back(node->create_subscription<MarkerArray>(
    raystarEndpoint(spawned.node_namespace, "non_homotopic_paths"), marker_qos,
    [&](MarkerArray::ConstSharedPtr message) {
      path_messages.push_back(std::move(message));
    }));
  subscriptions.push_back(node->create_subscription<MarkerArray>(
    raystarEndpoint(spawned.node_namespace, "poly_obstacles"), marker_qos,
    [&](MarkerArray::ConstSharedPtr message) {
      obstacle_messages.push_back(std::move(message));
    }));
  subscriptions.push_back(node->create_subscription<MarkerArray>(
    raystarEndpoint(spawned.node_namespace, "cdt"), marker_qos,
    [&](MarkerArray::ConstSharedPtr message) {
      cdt_messages.push_back(std::move(message));
    }));
  subscriptions.push_back(node->create_subscription<MarkerArray>(
    raystarEndpoint(spawned.node_namespace, "debug_tree"), marker_qos,
    [&](MarkerArray::ConstSharedPtr message) {
      tree_messages.push_back(std::move(message));
    }));

  const auto marker_stamp_ns = [](const MarkerArray& array) {
    const auto& stamp = array.markers.front().header.stamp;
    return static_cast<int64_t>(stamp.sec) * 1000000000LL +
      static_cast<int64_t>(stamp.nanosec);
  };
  const auto is_complete_replacement = [](const MarkerArray& array) {
    if (array.markers.size() < 2 ||
        array.markers.front().action !=
          visualization_msgs::msg::Marker::DELETEALL)
    {
      return false;
    }
    return std::all_of(array.markers.begin() + 1, array.markers.end(),
      [](const auto& marker) {
        return marker.action == visualization_msgs::msg::Marker::ADD;
      });
  };
  const auto is_replacement_snapshot = [&](const MarkerArray& array,
      int64_t min_stamp_ns, int64_t max_stamp_ns) {
    if (!is_complete_replacement(array))
      return false;
    const int64_t stamp_ns = marker_stamp_ns(array);
    if (stamp_ns < min_stamp_ns || stamp_ns > max_stamp_ns)
      return false;
    return true;
  };
  const auto is_clear_snapshot = [&](const MarkerArray& array,
      int64_t min_stamp_ns, int64_t max_stamp_ns) {
    return array.markers.size() == 1 &&
      array.markers.front().action ==
        visualization_msgs::msg::Marker::DELETEALL &&
      marker_stamp_ns(array) >= min_stamp_ns &&
      marker_stamp_ns(array) <= max_stamp_ns;
  };
  const auto find_message = [](const MarkerHistory& history, const auto& predicate)
    -> MarkerArray::ConstSharedPtr
  {
    const auto found = std::find_if(history.begin(), history.end(),
      [&](const auto& message) { return predicate(*message); });
    return found == history.end() ? nullptr : *found;
  };
  const auto spin_until = [&](const auto& predicate,
      std::chrono::milliseconds timeout) {
    const auto end_time = std::chrono::steady_clock::now() + timeout;
    while (std::chrono::steady_clock::now() < end_time)
    {
      executor.spin_some(50ms);
      if (predicate())
        return true;
      std::this_thread::sleep_for(10ms);
    }
    executor.spin_some(50ms);
    return predicate();
  };

  MarkerArray::ConstSharedPtr retained_paths;
  MarkerArray::ConstSharedPtr retained_obstacles;
  MarkerArray::ConstSharedPtr retained_cdt;
  MarkerArray::ConstSharedPtr retained_tree;
  ASSERT_TRUE(spin_until([&]() {
    retained_paths = find_message(path_messages, [&](const auto& array) {
      return is_complete_replacement(array) &&
        array.markers.size() == first_response->path_results.size() + 1;
    });
    retained_obstacles = find_message(obstacle_messages, [&](const auto& array) {
      return is_replacement_snapshot(
        array, first_request_start_ns, subscription_start_ns);
    });
    retained_cdt = find_message(cdt_messages, [&](const auto& array) {
      return is_replacement_snapshot(
        array, first_request_start_ns, subscription_start_ns);
    });
    retained_tree = find_message(tree_messages, [&](const auto& array) {
      return is_replacement_snapshot(
        array, first_request_start_ns, subscription_start_ns);
    });
    return retained_paths && retained_obstacles && retained_cdt && retained_tree;
  }, 3s)) << "Late subscribers did not receive retained marker snapshots";

  for (const auto* snapshot : {
      retained_paths.get(), retained_obstacles.get(),
      retained_cdt.get(), retained_tree.get()})
  {
    ASSERT_NE(snapshot, nullptr);
    for (const auto& marker : snapshot->markers)
      EXPECT_EQ(marker.header.frame_id, "map");
  }

  // RViz deletes a disabled path namespace locally and needs a fresh message
  // before it can display that namespace again.  Only the already-built path
  // snapshot should be repeated; the much larger geometry/debug snapshots
  // must remain one-shot transient-local publications.
  path_messages.clear();
  obstacle_messages.clear();
  cdt_messages.clear();
  tree_messages.clear();
  ASSERT_TRUE(spin_until([&]() {
    return path_messages.size() >= 2;
  }, 1500ms)) << "Cached path visualization was not periodically republished";
  EXPECT_TRUE(obstacle_messages.empty());
  EXPECT_TRUE(cdt_messages.empty());
  EXPECT_TRUE(tree_messages.empty());
  for (const auto& repeated_paths : path_messages)
  {
    ASSERT_NE(repeated_paths, nullptr);
    EXPECT_TRUE(*repeated_paths == *retained_paths)
      << "Periodic path publication was rebuilt instead of reusing the cache";
  }

  path_messages.clear();
  obstacle_messages.clear();
  cdt_messages.clear();
  tree_messages.clear();
  auto fewer_request = makeTestRequest();
  std::fill(fewer_request->map.data.begin(), fewer_request->map.data.end(), 0);
  fewer_request->k = 1;
  auto fewer_response = callService(executor, client, fewer_request);
  ASSERT_NE(fewer_response, nullptr);
  ASSERT_TRUE(fewer_response->success) << fewer_response->message;
  ASSERT_EQ(fewer_response->path_results.size(), 1u);

  MarkerArray::ConstSharedPtr fewer_paths;
  ASSERT_TRUE(spin_until([&]() {
    fewer_paths = find_message(path_messages, [&](const auto& array) {
      return is_complete_replacement(array) &&
        array.markers.size() == 2;
    });
    return fewer_paths != nullptr;
  }, 2s));
  ASSERT_EQ(fewer_paths->markers.size(), 2u);
  EXPECT_EQ(fewer_paths->markers[1].ns, "path_1");
  EXPECT_EQ(fewer_paths->markers[1].id, 1);

  // The periodic cache must be replaced transactionally with the new K=1
  // result rather than continuing to replay the previous K=3 snapshot.
  path_messages.clear();
  ASSERT_TRUE(spin_until([&]() {
    return !path_messages.empty();
  }, 1s));
  for (const auto& repeated_paths : path_messages)
  {
    ASSERT_NE(repeated_paths, nullptr);
    EXPECT_TRUE(*repeated_paths == *fewer_paths);
  }

  path_messages.clear();
  obstacle_messages.clear();
  cdt_messages.clear();
  tree_messages.clear();
  const int64_t invalid_request_start_ns = node->now().nanoseconds();
  auto invalid_request = makeTestRequest();
  invalid_request->map.data.pop_back();
  auto invalid_response = callService(executor, client, invalid_request);
  ASSERT_NE(invalid_response, nullptr);
  EXPECT_FALSE(invalid_response->success);
  const int64_t invalid_response_end_ns = node->now().nanoseconds();

  const auto has_clear_for_request = [&](const MarkerHistory& history) {
    return find_message(history, [&](const auto& array) {
      return is_clear_snapshot(
        array, invalid_request_start_ns, invalid_response_end_ns);
    }) != nullptr;
  };
  ASSERT_TRUE(spin_until([&]() {
    return has_clear_for_request(path_messages) &&
      has_clear_for_request(obstacle_messages) &&
      has_clear_for_request(cdt_messages) &&
      has_clear_for_request(tree_messages);
  }, 2s));

  // Observe by receive order, not by Marker header stamp: cached snapshots
  // deliberately keep their original stamp, so a stale replay after this
  // DELETEALL would otherwise evade a timestamp-based assertion.
  path_messages.clear();
  obstacle_messages.clear();
  cdt_messages.clear();
  tree_messages.clear();
  const auto no_resurrection_end = std::chrono::steady_clock::now() + 800ms;
  while (std::chrono::steady_clock::now() < no_resurrection_end)
  {
    executor.spin_some(50ms);
    std::this_thread::sleep_for(10ms);
  }
  const auto contains_add_after_failure = [](const MarkerHistory& history) {
    return std::any_of(history.begin(), history.end(), [](const auto& message) {
      return std::any_of(message->markers.begin(), message->markers.end(),
        [](const auto& marker) {
          return marker.action == visualization_msgs::msg::Marker::ADD;
        });
    });
  };
  EXPECT_FALSE(contains_add_after_failure(path_messages));
  EXPECT_FALSE(contains_add_after_failure(obstacle_messages));
  EXPECT_FALSE(contains_add_after_failure(cdt_messages));
  EXPECT_FALSE(contains_add_after_failure(tree_messages));
}

TEST_F(IntegrationTestFixture, ZeroPathRepublishPeriodKeepsOnlyDurableSnapshot)
{
  const auto spawned = spawnRaystarNode(
    "marker_refresh_disabled", {"path_visualization_republish_period_ms:=0"});
  ASSERT_GT(spawned.pid, 0)
    << "Unable to exec marker-refresh raystar_node: "
    << (spawned.launch_errno ? std::strerror(spawned.launch_errno) :
    "fork/pipe failure");
  ChildProcessGuard marker_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>(
    "test_marker_refresh_disabled_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);
  auto response = callService(executor, client, makeTestRequest());
  ASSERT_NE(response, nullptr);
  ASSERT_TRUE(response->success) << response->message;

  using MarkerArray = visualization_msgs::msg::MarkerArray;
  std::vector<MarkerArray::ConstSharedPtr> path_messages;
  auto subscription = node->create_subscription<MarkerArray>(
    raystarEndpoint(spawned.node_namespace, "non_homotopic_paths"),
    rclcpp::QoS(rclcpp::KeepLast(10)).transient_local(),
    [&](MarkerArray::ConstSharedPtr message) {
      path_messages.push_back(std::move(message));
    });
  (void)subscription;

  const auto retained_deadline = std::chrono::steady_clock::now() + 3s;
  while (path_messages.empty() &&
      std::chrono::steady_clock::now() < retained_deadline)
  {
    executor.spin_some(50ms);
    std::this_thread::sleep_for(10ms);
  }
  ASSERT_FALSE(path_messages.empty())
    << "Late subscriber did not receive a retained path snapshot";
  ASSERT_GT(path_messages.back()->markers.size(), 1u);
  EXPECT_EQ(path_messages.back()->markers.front().action,
    visualization_msgs::msg::Marker::DELETEALL);

  // Allow any discovery-time retained delivery already in flight to drain,
  // then require the observed count to stay fixed for longer than the normal
  // two-second default period.  This is robust to an RMW delivering the same
  // retained sample more than once during endpoint matching.
  const auto drain_deadline = std::chrono::steady_clock::now() + 200ms;
  while (std::chrono::steady_clock::now() < drain_deadline)
  {
    executor.spin_some(50ms);
    std::this_thread::sleep_for(10ms);
  }
  const size_t retained_count = path_messages.size();
  const auto no_refresh_deadline = std::chrono::steady_clock::now() + 2300ms;
  while (std::chrono::steady_clock::now() < no_refresh_deadline)
  {
    executor.spin_some(50ms);
    std::this_thread::sleep_for(10ms);
  }
  EXPECT_EQ(path_messages.size(), retained_count)
    << "A zero republish period must disable periodic path traffic";
}

TEST_F(IntegrationTestFixture, PathRepublishPeriodDefaultsToTwoSeconds)
{
  auto node = std::make_shared<rclcpp::Node>(
    "test_path_refresh_default_parameter_client");
  auto client = std::make_shared<rclcpp::AsyncParametersClient>(
    node, raystarEndpoint(main_namespace_));
  ASSERT_TRUE(client->wait_for_service(5s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);
  auto future = client->get_parameters(
    {"path_visualization_republish_period_ms"});
  ASSERT_TRUE(waitForFuture(executor, future, 5s));
  const auto parameters = future.get();
  ASSERT_EQ(parameters.size(), 1u);
  EXPECT_EQ(parameters.front().get_type(),
    rclcpp::ParameterType::PARAMETER_INTEGER);
  EXPECT_EQ(parameters.front().as_int(), 2000);
}

TEST_F(IntegrationTestFixture, CustomMapFramePropagatesToEveryReturnedPathPose)
{
  auto node = std::make_shared<rclcpp::Node>("test_custom_frame_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  constexpr const char* custom_frame = "warehouse_map";
  auto request = makeTestRequest();
  request->map.header.frame_id = custom_frame;
  request->start.header.frame_id = custom_frame;
  request->goal.header.frame_id = custom_frame;

  auto response = callService(executor, client, request);
  ASSERT_NE(response, nullptr);
  ASSERT_TRUE(response->success) << response->message;
  ASSERT_FALSE(response->path_results.empty());
  for (const auto& path_result : response->path_results)
  {
    const auto& path = path_result.path;
    EXPECT_EQ(path.header.frame_id, custom_frame);
    ASSERT_FALSE(path.poses.empty());
    for (const auto& pose : path.poses)
    {
      EXPECT_EQ(pose.header.frame_id, custom_frame);
      EXPECT_TRUE(std::isfinite(pose.pose.position.x));
      EXPECT_TRUE(std::isfinite(pose.pose.position.y));
    }
  }
}

TEST_F(IntegrationTestFixture, NegativeIdentityMapQuaternionIsAccepted)
{
  auto node = std::make_shared<rclcpp::Node>("test_equivalent_identity_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto request = makeTestRequest();
  request->map.info.origin.orientation.w = -1.0;
  auto response = callService(executor, client, request);
  ASSERT_NE(response, nullptr);
  EXPECT_TRUE(response->success) << response->message;
  EXPECT_FALSE(response->path_results.empty());
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
