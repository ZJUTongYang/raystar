#include <gtest/gtest.h>
#include <rclcpp/rclcpp.hpp>
#include <rclcpp_action/rclcpp_action.hpp>
#include <rclcpp/parameter_client.hpp>
#include <rcl_interfaces/msg/parameter_descriptor.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <visualization_msgs/msg/marker_array.hpp>
#include <raystar_interfaces/action/build_raystar_transition_graph.hpp>
#include <raystar_interfaces/action/plan_raystar_goal_set.hpp>
#include <raystar_interfaces/action/plan_raystar_paths.hpp>
#include <raystar_interfaces/environment_identity.hpp>
#include <raystar_interfaces/map_identity.hpp>
#include <raystar_interfaces/msg/map_status.hpp>
#include <raystar_interfaces/msg/planning_result_info.hpp>
#include <raystar_interfaces/srv/get_raystar_paths.hpp>

#include "path_collision_oracle.h"
#include "small_grid_property.h"
#include "conservative_path_length.h"

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
using RaystarGoalSetAction = raystar_interfaces::action::PlanRaystarGoalSet;
using RaystarTransitionAction = raystar_interfaces::action::BuildRaystarTransitionGraph;
using PlanningResultInfo = raystar_interfaces::msg::PlanningResultInfo;

static nav_msgs::msg::OccupancyGrid makeTestGrid() {
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
    for (unsigned int x = 10; x < 20; ++x) grid.data[y * 30 + x] = 100;

  return grid;
}

static nav_msgs::msg::OccupancyGrid makeRawContourWedgeMap() {
  nav_msgs::msg::OccupancyGrid grid;
  grid.header.frame_id = "map";
  grid.info.width = 50;
  grid.info.height = 50;
  grid.info.resolution = 1.0f;
  grid.info.origin.orientation.w = 1.0;
  grid.data.assign(50U * 50U, 0);

  for (unsigned int x = 0; x < grid.info.width; ++x) {
    grid.data[x] = 100;
    grid.data[(grid.info.height - 1U) * grid.info.width + x] = 100;
  }
  for (unsigned int y = 0; y < grid.info.height; ++y) {
    grid.data[y * grid.info.width] = 100;
    grid.data[y * grid.info.width + grid.info.width - 1U] = 100;
  }
  return grid;
}

static RaystarService::Request::SharedPtr makeTestRequest() {
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

static RaystarService::Request::SharedPtr makeVerticalBarrierRequest(int8_t occupancy_value,
                                                                     bool allow_unknown = false) {
  auto request = makeTestRequest();
  std::fill(request->map.data.begin(), request->map.data.end(), 0);
  const size_t width = static_cast<size_t>(request->map.info.width);
  const size_t height = static_cast<size_t>(request->map.info.height);
  const size_t barrier_x = width / 2;
  for (size_t y = 0; y < height; ++y) request->map.data[y * width + barrier_x] = occupancy_value;
  request->k = 1;
  request->allow_unknown = allow_unknown;
  return request;
}

static RaystarAction::Goal makeTestActionGoal() {
  auto request = makeTestRequest();
  RaystarAction::Goal goal;
  goal.map_id = raystar_interfaces::computeMapId(request->map);
  goal.start = std::move(request->start);
  goal.goal = std::move(request->goal);
  goal.search_mode = request->search_mode;
  goal.k = request->k;
  goal.max_path_length = request->max_path_length;
  goal.allow_self_crossing = request->allow_self_crossing;
  goal.allow_unknown = request->allow_unknown;
  goal.include_debug = false;
  return goal;
}

static RaystarGoalSetAction::Goal makeTestGoalSetActionGoal() {
  auto request = makeTestRequest();
  RaystarGoalSetAction::Goal goal;
  goal.map_id = raystar_interfaces::computeMapId(request->map);
  goal.start = request->start;
  auto upper = request->goal;
  upper.pose.position.y = 14.5;
  auto lower = request->goal;
  lower.pose.position.y = 16.5;
  goal.goals = {upper, lower};
  goal.max_path_lengths = {45.0, 45.0};
  goal.allow_self_crossing = false;
  goal.allow_unknown = false;
  goal.include_debug = false;
  return goal;
}

static nav_msgs::msg::OccupancyGrid makeLongRunningActionMap() {
  constexpr unsigned int width = 768;
  constexpr unsigned int height = 768;

  auto map = makeTestGrid();
  map.info.width = width;
  map.info.height = height;
  map.data.assign(static_cast<size_t>(width) * height, 0);

  // Thousands of isolated one-cell obstacles make contour/CDT construction
  // decisively longer than the Action goal/cancel round trips, while keeping
  // one connected free-space component and avoiding diagonal vertex contact.
  for (unsigned int y = 4; y + 4 < height; y += 4) {
    for (unsigned int x = 4; x + 4 < width; x += 4)
      map.data[static_cast<size_t>(y) * width + x] = 100;
  }
  return map;
}

static RaystarAction::Goal makeLongRunningActionGoal() {
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

template <typename FutureT>
static bool waitForFuture(rclcpp::executors::SingleThreadedExecutor& executor,
                          FutureT& future,
                          std::chrono::milliseconds timeout) {
  const auto deadline = std::chrono::steady_clock::now() + timeout;
  while (std::chrono::steady_clock::now() < deadline) {
    if (future.wait_for(0s) == std::future_status::ready)
      return true;
    executor.spin_some(20ms);
  }
  executor.spin_some(20ms);
  return future.wait_for(0s) == std::future_status::ready;
}

static bool cacheMapAndWait(rclcpp::executors::SingleThreadedExecutor& executor,
                            const rclcpp::Node::SharedPtr& node,
                            const nav_msgs::msg::OccupancyGrid& map,
                            const std::string& server_namespace,
                            std::chrono::milliseconds timeout = 5s) {
  const auto expected_id = raystar_interfaces::computeMapId(map);
  const auto expected_environment_id_disallow_unknown =
    raystar_interfaces::computeEnvironmentId(map, 99, false);
  const auto expected_environment_id_allow_unknown =
    raystar_interfaces::computeEnvironmentId(map, 99, true);
  const auto status_topic = server_namespace + "/raystar/map_status";
  const auto map_topic = server_namespace + "/map";
  bool admitted = false;
  bool saw_matching_map_status = false;
  bool environment_ids_nonzero = false;
  auto status_subscription = node->create_subscription<raystar_interfaces::msg::MapStatus>(
    status_topic,
    rclcpp::QoS(rclcpp::KeepLast(1)).transient_local(),
    [&](raystar_interfaces::msg::MapStatus::ConstSharedPtr status) {
      if (!status || !raystar_interfaces::mapIdsEqual(status->map_id, expected_id)) {
        return;
      }
      saw_matching_map_status = true;
      environment_ids_nonzero =
        !raystar_interfaces::isZeroEnvironmentId(status->environment_id_disallow_unknown) &&
        !raystar_interfaces::isZeroEnvironmentId(status->environment_id_allow_unknown);
      admitted =
        environment_ids_nonzero &&
        status->occupied_threshold == 99U &&
        status->environment_identity_version == raystar_interfaces::kEnvironmentIdentityVersion &&
        status->occupancy_semantics_version == raystar_interfaces::kOccupancySemanticsVersion &&
        status->geometry_semantics_version == raystar_interfaces::kGeometrySemanticsVersion &&
        status->topology_semantics_version == raystar_interfaces::kTopologySemanticsVersion &&
        raystar_interfaces::environmentIdsEqual(status->environment_id_disallow_unknown,
                                                expected_environment_id_disallow_unknown) &&
        raystar_interfaces::environmentIdsEqual(status->environment_id_allow_unknown,
                                                expected_environment_id_allow_unknown);
    });
  auto map_publisher = node->create_publisher<nav_msgs::msg::OccupancyGrid>(
    map_topic, rclcpp::QoS(rclcpp::KeepLast(1)).transient_local());

  const auto deadline = std::chrono::steady_clock::now() + timeout;
  auto next_publish = std::chrono::steady_clock::time_point::min();
  while (!admitted && std::chrono::steady_clock::now() < deadline) {
    const auto now = std::chrono::steady_clock::now();
    if (now >= next_publish) {
      map_publisher->publish(map);
      next_publish = now + 100ms;
    }
    executor.spin_some(20ms);
    std::this_thread::sleep_for(1ms);
  }
  (void)status_subscription;
  if (saw_matching_map_status) {
    EXPECT_TRUE(environment_ids_nonzero)
      << "MapStatus environment identities must both be nonzero";
  }
  return admitted;
}

// Every integration-test invocation gets a private ROS namespace.  This is
// intentionally based on the parent PID rather than a fixed literal so that
// two ctest jobs can run concurrently without sharing services, actions or
// cached-map topics.  Child nodes inherit the same namespace prefix.
static const std::string& integrationNamespacePrefix() {
  static const std::string prefix =
    "/raystar_integration_" + std::to_string(static_cast<long long>(getpid()));
  return prefix;
}

static std::string integrationNamespace(const std::string& scenario) {
  return integrationNamespacePrefix() + "/" + scenario;
}

static std::string raystarEndpoint(const std::string& server_namespace,
                                   const std::string& suffix = {}) {
  return server_namespace + "/raystar" + (suffix.empty() ? std::string{} : "/" + suffix);
}

static std::string mainServerNamespace() {
  return integrationNamespace("main");
}

struct SpawnedNode {
  pid_t pid{-1};
  std::string node_namespace;
  int launch_errno{0};
};

// Start raystar_node and synchronously distinguish a successful exec from an
// immediate child-side failure.  The close-on-exec pipe reaches EOF only after
// exec has replaced the child image; otherwise the child writes errno before
// exiting.  This avoids reporting a discovery timeout when the executable (or
// its dynamic loader) was actually missing.
static SpawnedNode spawnRaystarNode(const std::string& scenario,
                                    const std::vector<std::string>& parameter_overrides = {}) {
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
  for (const auto& override : parameter_overrides) {
    arguments.emplace_back("-p");
    arguments.push_back(override);
  }

  int error_pipe[2] = {-1, -1};
  if (pipe(error_pipe) != 0) {
    spawned.launch_errno = errno;
    return spawned;
  }
  const int flags = fcntl(error_pipe[1], F_GETFD);
  if (flags < 0 || fcntl(error_pipe[1], F_SETFD, flags | FD_CLOEXEC) < 0) {
    spawned.launch_errno = errno;
    close(error_pipe[0]);
    close(error_pipe[1]);
    return spawned;
  }

  std::vector<char*> argv;
  argv.reserve(arguments.size() + 1);
  for (auto& argument : arguments) argv.push_back(argument.data());
  argv.push_back(nullptr);

  spawned.pid = fork();
  if (spawned.pid < 0) {
    spawned.launch_errno = errno;
    close(error_pipe[0]);
    close(error_pipe[1]);
    return spawned;
  }
  if (spawned.pid == 0) {
    close(error_pipe[0]);
    execv(argv[0], argv.data());
    const int child_errno = errno;
    // write(2) is async-signal-safe and therefore suitable in the post-fork
    // child even though the parent process owns DDS/background threads.
    const char* bytes = reinterpret_cast<const char*>(&child_errno);
    size_t written = 0;
    while (written < sizeof(child_errno)) {
      const ssize_t count = write(error_pipe[1], bytes + written, sizeof(child_errno) - written);
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
  while (received < sizeof(child_errno)) {
    const ssize_t count = read(error_pipe[0],
                               reinterpret_cast<char*>(&child_errno) + received,
                               sizeof(child_errno) - received);
    if (count > 0)
      received += static_cast<size_t>(count);
    else if (count < 0 && errno == EINTR)
      continue;
    else
      break;
  }
  close(error_pipe[0]);
  if (received != 0) {
    spawned.launch_errno = child_errno;
    int ignored_status = 0;
    (void)waitpid(spawned.pid, &ignored_status, 0);
    spawned.pid = -1;
  }
  return spawned;
}

static bool waitForChildExit(pid_t pid, int& status, std::chrono::milliseconds timeout) {
  const auto deadline = std::chrono::steady_clock::now() + timeout;
  while (std::chrono::steady_clock::now() < deadline) {
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
  const RaystarService::Request::SharedPtr& request) {
  auto future = client->async_send_request(request);
  const auto end_time = std::chrono::steady_clock::now() + 10s;
  while (std::chrono::steady_clock::now() < end_time) {
    if (future.wait_for(100ms) == std::future_status::ready)
      break;
    executor.spin_some(100ms);
  }
  if (future.wait_for(0s) != std::future_status::ready) {
    ADD_FAILURE() << "Service call did not complete within 10 seconds";
    return nullptr;
  }
  return future.get();
}

static std::string lowercase(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(), [](unsigned char character) {
    return static_cast<char>(std::tolower(character));
  });
  return text;
}

template <typename ResponseT>
static void expectStructuredResultInvariants(const ResponseT& response) {
  const auto& info = response.result_info;
  EXPECT_EQ(info.returned_path_count, response.path_results.size());
  EXPECT_LE(info.returned_path_count, info.found_path_count);
  if (info.search_mode == RaystarAction::Goal::SEARCH_MODE_TOP_K)
    EXPECT_LE(info.found_path_count, info.requested_path_count);
  else
    EXPECT_EQ(info.requested_path_count, 0u);
  EXPECT_EQ(response.success, !response.path_results.empty());
  if (info.request_satisfied) {
    EXPECT_TRUE(info.search_complete);
    EXPECT_TRUE(info.output_complete);
    if (info.search_mode == RaystarAction::Goal::SEARCH_MODE_TOP_K)
      EXPECT_EQ(info.returned_path_count, info.requested_path_count);
    else
      EXPECT_TRUE(info.cost_bound_exhausted);
  }
}

static void expectIndependentCollisionFreePath(const nav_msgs::msg::OccupancyGrid& map,
                                               const nav_msgs::msg::Path& path,
                                               bool allow_unknown,
                                               int occupied_threshold = 99) {
  raystar::test_oracle::OracleOptions options;
  options.occupied_threshold = occupied_threshold;
  options.allow_unknown = allow_unknown;
  const auto result = raystar::test_oracle::validatePath(map, path, options);
  EXPECT_EQ(result.status, raystar::test_oracle::ValidationStatus::kCollisionFree)
    << result.diagnostic;
}

static void expectIndependentSelfIntersectionFreePath(const nav_msgs::msg::Path& path) {
  const auto result = raystar::test_oracle::validateNoSelfIntersection(path);
  EXPECT_EQ(result.status, raystar::test_oracle::SelfIntersectionStatus::kIntersectionFree)
    << result.diagnostic;
}

static RaystarService::Request::SharedPtr makeSmallGridPropertyRequest(
  const raystar::test_property::SmallGridCase& property_case) {
  auto request = std::make_shared<RaystarService::Request>();
  request->map = property_case.map;
  request->start.header.frame_id = property_case.map.header.frame_id;
  request->goal.header.frame_id = property_case.map.header.frame_id;

  const double resolution = static_cast<double>(property_case.map.info.resolution);
  const auto cell_center = [&](const raystar::test_property::GridCell& cell) {
    return std::pair<double, double>{
      property_case.map.info.origin.position.x + (static_cast<double>(cell.x) + 0.5) * resolution,
      property_case.map.info.origin.position.y + (static_cast<double>(cell.y) + 0.5) * resolution};
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

static double polylineMetricLength(const nav_msgs::msg::Path& path) {
  double length = 0.0;
  for (size_t index = 1; index < path.poses.size(); ++index) {
    const auto& previous = path.poses[index - 1].pose.position;
    const auto& current = path.poses[index].pose.position;
    length += std::hypot(current.x - previous.x, current.y - previous.y);
  }
  return length;
}

static void expectSmallGridPathContract(const raystar::test_property::SmallGridCase& property_case,
                                        const RaystarService::Request& request,
                                        const raystar_interfaces::msg::PathResult& path_result) {
  const auto& path = path_result.path;
  const auto& topology_path = path_result.topology_path;
  EXPECT_EQ(path.header.frame_id, property_case.map.header.frame_id);
  EXPECT_EQ(topology_path.header.frame_id, property_case.map.header.frame_id);
  ASSERT_GE(path.poses.size(), 2u);
  ASSERT_GE(topology_path.poses.size(), 2u);
  for (size_t pose_index = 0; pose_index < path.poses.size(); ++pose_index) {
    SCOPED_TRACE("pose index " + std::to_string(pose_index));
    const auto& pose_stamped = path.poses[pose_index];
    const auto& position = pose_stamped.pose.position;
    const auto& orientation = pose_stamped.pose.orientation;
    EXPECT_EQ(pose_stamped.header.frame_id, property_case.map.header.frame_id);
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
  const auto& topology_first = topology_path.poses.front().pose.position;
  const auto& topology_last = topology_path.poses.back().pose.position;
  EXPECT_NEAR(first.x, request.start.pose.position.x, 1e-9);
  EXPECT_NEAR(first.y, request.start.pose.position.y, 1e-9);
  EXPECT_NEAR(last.x, request.goal.pose.position.x, 1e-9);
  EXPECT_NEAR(last.y, request.goal.pose.position.y, 1e-9);
  EXPECT_NEAR(topology_first.x, request.start.pose.position.x, 1e-9);
  EXPECT_NEAR(topology_first.y, request.start.pose.position.y, 1e-9);
  EXPECT_NEAR(topology_last.x, request.goal.pose.position.x, 1e-9);
  EXPECT_NEAR(topology_last.y, request.goal.pose.position.y, 1e-9);

  EXPECT_TRUE(std::isfinite(path_result.cost));
  EXPECT_GE(path_result.cost, 0.0);
  const double sampled_length = polylineMetricLength(path);
  const double topology_length = polylineMetricLength(topology_path);
  const double length_tolerance = 1e-8 * std::max(1.0, sampled_length);
  EXPECT_NEAR(path_result.cost, sampled_length, length_tolerance);
  EXPECT_NEAR(path_result.cost, topology_length, length_tolerance);
  const double direct_distance =
    std::hypot(request.goal.pose.position.x - request.start.pose.position.x,
               request.goal.pose.position.y - request.start.pose.position.y);
  EXPECT_GE(path_result.cost + length_tolerance, direct_distance);
  expectIndependentCollisionFreePath(property_case.map,
                                     path,
                                     property_case.allow_unknown,
                                     raystar::test_property::kOccupiedThreshold);
  expectIndependentCollisionFreePath(property_case.map,
                                     topology_path,
                                     property_case.allow_unknown,
                                     raystar::test_property::kOccupiedThreshold);
  expectIndependentSelfIntersectionFreePath(path);
  expectIndependentSelfIntersectionFreePath(topology_path);
}

static void expectSmallGridReachableResult(
  const raystar::test_property::SmallGridCase& property_case,
  const RaystarService::Request& request,
  const RaystarService::Response& response) {
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
  EXPECT_TRUE(raystar_interfaces::mapIdsEqual(info.map_id,
                                              raystar_interfaces::computeMapId(property_case.map)));
  EXPECT_TRUE(raystar_interfaces::environmentIdsEqual(
    info.environment_id,
    raystar_interfaces::computeEnvironmentId(
      property_case.map, raystar::test_property::kOccupiedThreshold, property_case.allow_unknown)));
  EXPECT_TRUE(std::isfinite(info.map_time_ms));
  EXPECT_GE(info.map_time_ms, 0.0);
  EXPECT_TRUE(std::isfinite(info.plan_time_ms));
  EXPECT_GE(info.plan_time_ms, 0.0);
  EXPECT_TRUE(response.debug_nodes.empty());
  ASSERT_EQ(response.path_results.size(), 1u);
  expectSmallGridPathContract(property_case, request, response.path_results.front());
  expectStructuredResultInvariants(response);
}

static void expectSmallGridUnreachableResult(
  const raystar::test_property::SmallGridCase& property_case,
  const RaystarService::Response& response) {
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
  EXPECT_TRUE(raystar_interfaces::mapIdsEqual(info.map_id,
                                              raystar_interfaces::computeMapId(property_case.map)));
  EXPECT_TRUE(raystar_interfaces::environmentIdsEqual(
    info.environment_id,
    raystar_interfaces::computeEnvironmentId(
      property_case.map, raystar::test_property::kOccupiedThreshold, property_case.allow_unknown)));
  EXPECT_TRUE(std::isfinite(info.map_time_ms));
  EXPECT_GE(info.map_time_ms, 0.0);
  EXPECT_TRUE(std::isfinite(info.plan_time_ms));
  EXPECT_GE(info.plan_time_ms, 0.0);
  EXPECT_TRUE(response.path_results.empty());
  EXPECT_TRUE(response.debug_nodes.empty());
  expectStructuredResultInvariants(response);
}

static void expectStableSmallGridResponses(const RaystarService::Response& first,
                                           const RaystarService::Response& second) {
  EXPECT_EQ(first.success, second.success);
  const auto& first_info = first.result_info;
  const auto& second_info = second.result_info;
  EXPECT_TRUE(raystar_interfaces::mapIdsEqual(first_info.map_id, second_info.map_id));
  EXPECT_TRUE(
    raystar_interfaces::environmentIdsEqual(first_info.environment_id, second_info.environment_id));
  EXPECT_EQ(first_info.status, second_info.status);
  EXPECT_EQ(first_info.limits_reached, second_info.limits_reached);
  EXPECT_EQ(first_info.request_satisfied, second_info.request_satisfied);
  EXPECT_EQ(first_info.search_complete, second_info.search_complete);
  EXPECT_EQ(first_info.output_complete, second_info.output_complete);
  EXPECT_EQ(first_info.debug_requested, second_info.debug_requested);
  EXPECT_EQ(first_info.debug_output_complete, second_info.debug_output_complete);
  EXPECT_EQ(first_info.requested_path_count, second_info.requested_path_count);
  EXPECT_EQ(first_info.found_path_count, second_info.found_path_count);
  EXPECT_EQ(first_info.returned_path_count, second_info.returned_path_count);
  EXPECT_EQ(first_info.search_mode, second_info.search_mode);
  EXPECT_DOUBLE_EQ(first_info.requested_max_path_length, second_info.requested_max_path_length);
  EXPECT_EQ(first_info.cost_bound_exhausted, second_info.cost_bound_exhausted);
  EXPECT_EQ(first.path_results.size(), second.path_results.size());
  EXPECT_EQ(first.debug_nodes.size(), second.debug_nodes.size());

  ASSERT_EQ(first.path_results.size(), second.path_results.size());
  for (size_t path_index = 0; path_index < first.path_results.size(); ++path_index) {
    SCOPED_TRACE("repeated path index " + std::to_string(path_index));
    const auto& first_result = first.path_results[path_index];
    const auto& second_result = second.path_results[path_index];
    EXPECT_DOUBLE_EQ(first_result.cost, second_result.cost);
    EXPECT_EQ(first_result.path.header.frame_id, second_result.path.header.frame_id);
    ASSERT_EQ(first_result.path.poses.size(), second_result.path.poses.size());
    for (size_t pose_index = 0; pose_index < first_result.path.poses.size(); ++pose_index) {
      SCOPED_TRACE("repeated pose index " + std::to_string(pose_index));
      const auto& first_pose = first_result.path.poses[pose_index];
      const auto& second_pose = second_result.path.poses[pose_index];
      EXPECT_EQ(first_pose.header.frame_id, second_pose.header.frame_id);
      EXPECT_DOUBLE_EQ(first_pose.pose.position.x, second_pose.pose.position.x);
      EXPECT_DOUBLE_EQ(first_pose.pose.position.y, second_pose.pose.position.y);
      EXPECT_DOUBLE_EQ(first_pose.pose.position.z, second_pose.pose.position.z);
      EXPECT_DOUBLE_EQ(first_pose.pose.orientation.x, second_pose.pose.orientation.x);
      EXPECT_DOUBLE_EQ(first_pose.pose.orientation.y, second_pose.pose.orientation.y);
      EXPECT_DOUBLE_EQ(first_pose.pose.orientation.z, second_pose.pose.orientation.z);
      EXPECT_DOUBLE_EQ(first_pose.pose.orientation.w, second_pose.pose.orientation.w);
    }
  }
}

template <typename ResponseT>
static void expectSingleObstacleSimplePathContract(const nav_msgs::msg::OccupancyGrid& map,
                                                   const ResponseT& response) {
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
  EXPECT_TRUE(raystar_interfaces::mapIdsEqual(info.map_id, raystar_interfaces::computeMapId(map)));
  EXPECT_TRUE(raystar_interfaces::environmentIdsEqual(
    info.environment_id, raystar_interfaces::computeEnvironmentId(map, 99, false)));
  EXPECT_TRUE(std::isfinite(info.map_time_ms));
  EXPECT_GE(info.map_time_ms, 0.0);
  EXPECT_TRUE(std::isfinite(info.plan_time_ms));
  EXPECT_GE(info.plan_time_ms, 0.0);
  EXPECT_GT(info.expanded_nodes, 0u);
  EXPECT_EQ(response.path_results.size(), 2u);
  expectStructuredResultInvariants(response);

  for (size_t index = 0; index < response.path_results.size(); ++index) {
    SCOPED_TRACE("path index " + std::to_string(index));
    const auto& path = response.path_results[index].path;
    expectIndependentCollisionFreePath(map, path, false);
    expectIndependentSelfIntersectionFreePath(path);
  }
}

class ChildProcessGuard {
public:
  explicit ChildProcessGuard(pid_t pid) : pid_(pid) {}
  ChildProcessGuard(const ChildProcessGuard&) = delete;
  ChildProcessGuard& operator=(const ChildProcessGuard&) = delete;

  ~ChildProcessGuard() {
    stopAndValidate();
  }

  void disarm() noexcept {
    pid_ = -1;
  }

private:
  void stopAndValidate() noexcept {
    if (pid_ <= 0)
      return;

    int status = 0;
    bool reaped = waitForChildExit(pid_, status, 0ms);
    bool used_sigterm = false;
    bool used_sigkill = false;
    if (!reaped) {
      if (kill(pid_, SIGINT) < 0 && errno != ESRCH)
        ADD_FAILURE() << "Failed to send SIGINT to raystar_node: " << std::strerror(errno);
      reaped = waitForChildExit(pid_, status, 3s);
    }
    if (!reaped) {
      used_sigterm = true;
      if (kill(pid_, SIGTERM) < 0 && errno != ESRCH)
        ADD_FAILURE() << "Failed to send SIGTERM to raystar_node: " << std::strerror(errno);
      reaped = waitForChildExit(pid_, status, 2s);
    }
    if (!reaped) {
      used_sigkill = true;
      if (kill(pid_, SIGKILL) < 0 && errno != ESRCH)
        ADD_FAILURE() << "Failed to send SIGKILL to raystar_node: " << std::strerror(errno);
      reaped = waitForChildExit(pid_, status, 2s);
    }

    if (!reaped) {
      ADD_FAILURE() << "raystar_node child did not exit within bounded "
                    << "SIGINT/SIGTERM/SIGKILL cleanup";
      return;
    }
    if (used_sigkill)
      ADD_FAILURE() << "raystar_node required SIGKILL during cleanup";
    if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
      ADD_FAILURE() << "raystar_node exited unexpectedly (status=" << status << ")";
    }
    pid_ = -1;
    (void)used_sigterm;
  }

  pid_t pid_;
};

class IntegrationTestFixture : public ::testing::Test {
protected:
  static void SetUpTestSuite() {
    // Keep DDS discovery and rosout files isolated from a concurrently
    // running integration binary.  A caller-provided domain is retained so
    // CI can deliberately place tests in a preselected domain; otherwise a
    // PID-derived valid domain is assigned before rclcpp initializes.
    const auto existing_domain = std::getenv("ROS_DOMAIN_ID");
    const int domain_id =
      existing_domain ? std::atoi(existing_domain) : 20 + static_cast<int>(getpid() % 180);
    ASSERT_GE(domain_id, 0);
    ASSERT_LE(domain_id, 232);
    const std::string domain_text = std::to_string(domain_id);
    ASSERT_EQ(setenv("ROS_DOMAIN_ID", domain_text.c_str(), 1), 0) << std::strerror(errno);

    const std::filesystem::path log_directory =
      std::filesystem::temp_directory_path() /
      ("raystar_integration_" + std::to_string(static_cast<long long>(getpid())));
    std::error_code directory_error;
    std::filesystem::create_directories(log_directory, directory_error);
    ASSERT_FALSE(directory_error) << "Unable to create ROS_LOG_DIR " << log_directory << ": "
                                  << directory_error.message();
    ASSERT_EQ(setenv("ROS_LOG_DIR", log_directory.c_str(), 1), 0) << std::strerror(errno);

    rclcpp::init(0, nullptr);

    const auto spawned = spawnRaystarNode("main");
    ASSERT_GT(spawned.pid, 0) << "Unable to exec raystar_node: "
                              << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                       : "fork/pipe failure");
    pid_ = spawned.pid;
    main_namespace_ = spawned.node_namespace;

    auto node = std::make_shared<rclcpp::Node>("test_wait_node");
    auto client = node->create_client<raystar_interfaces::srv::GetRaystarPaths>(
      raystarEndpoint(main_namespace_, "get_raystar_paths"));
    auto action_client = rclcpp_action::create_client<RaystarAction>(
      node, raystarEndpoint(main_namespace_, "plan_paths"));

    bool found = false;
    for (int i = 0; i < 50; ++i) {
      if (client->wait_for_service(100ms)) {
        found = true;
        break;
      }
    }
    ASSERT_TRUE(found) << "Service " << raystarEndpoint(main_namespace_, "get_raystar_paths")
                       << " not available after 5s";
    ASSERT_TRUE(action_client->wait_for_action_server(5s))
      << "Action " << raystarEndpoint(main_namespace_, "plan_paths") << " not available after 5s";

    rclcpp::executors::SingleThreadedExecutor executor;
    executor.add_node(node);
    ASSERT_TRUE(cacheMapAndWait(executor, node, makeTestGrid(), main_namespace_))
      << "Raystar server did not admit the default cached map";
  }

  static void TearDownTestSuite() {
    if (pid_ > 0) {
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

TEST_F(IntegrationTestFixture, ServiceCallReturnsPaths) {
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
    if (status == std::future_status::ready)
      break;
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

TEST_F(IntegrationTestFixture, FinalPublicCostSortsAffineRoundingReversal) {
  auto node = std::make_shared<rclcpp::Node>("public_cost_order_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));

  auto request = std::make_shared<RaystarService::Request>();
  auto& map = request->map;
  map.header.frame_id = "map";
  map.info.width = 64;
  map.info.height = 64;
  map.info.resolution = 0.05f;
  const double resolution = static_cast<double>(map.info.resolution);
  // Shift the complete grid geometry one cell away from Raystar's reserved
  // outer ring, while shifting the origin back by that same binary64
  // resolution. The serialized world coordinates therefore retain the
  // audited affine-rounding reversal.
  map.info.origin.position.x = -resolution;
  map.info.origin.position.y = std::fma(-1.0, resolution, 0.1);
  map.info.origin.orientation.w = 1.0;
  map.data.assign(static_cast<size_t>(map.info.width) * map.info.height, 0);
  for (size_t y = 3; y < 14; ++y) {
    for (size_t x = 3; x < 14; ++x) map.data[y * map.info.width + x] = 100;
  }
  const auto world_x = [&](double grid_x) {
    return std::fma(grid_x, resolution, map.info.origin.position.x);
  };
  const auto world_y = [&](double grid_y) {
    return std::fma(grid_y, resolution, map.info.origin.position.y);
  };
  request->start.header.frame_id = "map";
  request->start.pose.position.x = world_x(1.5);
  request->start.pose.position.y = world_y(1.5);
  request->start.pose.orientation.w = 1.0;
  request->goal.header.frame_id = "map";
  request->goal.pose.position.x = world_x(32.5);
  request->goal.pose.position.y = world_y(32.5);
  request->goal.pose.orientation.w = 1.0;
  request->k = 2;

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);
  const auto response = callService(executor, client, request);

  ASSERT_NE(response, nullptr);
  ASSERT_TRUE(response->success) << response->message;
  ASSERT_EQ(response->path_results.size(), 2u) << response->message;
  EXPECT_DOUBLE_EQ(response->path_results[0].cost, 0x1.2f6d9bc5d1cdcp+1);
  EXPECT_DOUBLE_EQ(response->path_results[1].cost, 0x1.2f6d9bc5d1cddp+1);
  for (const auto& result : response->path_results)
    ASSERT_EQ(result.topology_path.poses.size(), 3u);
  const auto& lower_turn = response->path_results[0].topology_path.poses[1].pose.position;
  const auto& higher_turn = response->path_results[1].topology_path.poses[1].pose.position;
  EXPECT_DOUBLE_EQ(lower_turn.x, world_x(3.0));
  EXPECT_DOUBLE_EQ(lower_turn.y, world_y(14.0));
  EXPECT_DOUBLE_EQ(higher_turn.x, world_x(14.0));
  EXPECT_DOUBLE_EQ(higher_turn.y, world_y(3.0));
}

TEST_F(IntegrationTestFixture, GoalSetActionReturnsOrderedPerGoalCertificates) {
  auto node = std::make_shared<rclcpp::Node>("test_goal_set_action_client");
  auto client = rclcpp_action::create_client<RaystarGoalSetAction>(
    node, raystarEndpoint(main_namespace_, "plan_goal_set"));
  ASSERT_TRUE(client->wait_for_action_server(5s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto invalid_goal = makeTestGoalSetActionGoal();
  invalid_goal.max_path_lengths.pop_back();
  auto invalid_goal_future = client->async_send_goal(invalid_goal);
  ASSERT_TRUE(waitForFuture(executor, invalid_goal_future, 5s));
  const auto invalid_handle = invalid_goal_future.get();
  ASSERT_NE(invalid_handle, nullptr);
  auto invalid_result_future = client->async_get_result(invalid_handle);
  ASSERT_TRUE(waitForFuture(executor, invalid_result_future, 5s));
  const auto invalid_wrapped = invalid_result_future.get();
  EXPECT_EQ(invalid_wrapped.code, rclcpp_action::ResultCode::ABORTED);
  ASSERT_NE(invalid_wrapped.result, nullptr);
  EXPECT_EQ(invalid_wrapped.result->result_info.status, PlanningResultInfo::STATUS_INVALID_REQUEST);
  EXPECT_FALSE(invalid_wrapped.result->success);
  EXPECT_TRUE(invalid_wrapped.result->goal_results.empty());

  const auto valid_goal = makeTestGoalSetActionGoal();
  auto goal_future = client->async_send_goal(valid_goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto wrapped = result_future.get();

  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(wrapped.result, nullptr);
  EXPECT_TRUE(wrapped.result->success) << wrapped.result->message;
  EXPECT_EQ(wrapped.result->requested_start, valid_goal.start);
  EXPECT_EQ(wrapped.result->requested_allow_self_crossing, valid_goal.allow_self_crossing);
  EXPECT_EQ(wrapped.result->requested_allow_unknown, valid_goal.allow_unknown);
  EXPECT_EQ(wrapped.result->result_info.requested_goal_count, 2u);
  EXPECT_EQ(wrapped.result->result_info.returned_goal_count, 2u);
  EXPECT_EQ(wrapped.result->result_info.completed_goal_count, 2u);
  EXPECT_TRUE(wrapped.result->result_info.search_complete);
  EXPECT_TRUE(wrapped.result->result_info.output_complete);
  EXPECT_TRUE(wrapped.result->result_info.request_satisfied);
  const auto expected_environment_id =
    raystar_interfaces::computeEnvironmentId(makeTestGrid(), 99, false);
  EXPECT_TRUE(raystar_interfaces::environmentIdsEqual(wrapped.result->result_info.environment_id,
                                                      expected_environment_id));
  ASSERT_EQ(wrapped.result->goal_results.size(), 2u);
  for (size_t i = 0; i < wrapped.result->goal_results.size(); ++i) {
    const auto& goal_result = wrapped.result->goal_results[i];
    EXPECT_EQ(goal_result.goal_index, i);
    EXPECT_TRUE(goal_result.success) << goal_result.message;
    EXPECT_TRUE(goal_result.result_info.search_complete);
    EXPECT_TRUE(goal_result.result_info.cost_bound_exhausted);
    EXPECT_TRUE(goal_result.result_info.output_complete);
    EXPECT_TRUE(goal_result.result_info.request_satisfied);
    EXPECT_TRUE(raystar_interfaces::environmentIdsEqual(goal_result.result_info.environment_id,
                                                        expected_environment_id));
    EXPECT_FALSE(goal_result.path_results.empty());
    EXPECT_DOUBLE_EQ(goal_result.goal.pose.position.y, i == 0 ? 14.5 : 16.5);
  }
}

TEST_F(IntegrationTestFixture, GoalSetCompleteEmptyEnumerationUsesSucceededTransport) {
  auto node = std::make_shared<rclcpp::Node>("test_empty_goal_set_action_client");
  auto client = rclcpp_action::create_client<RaystarGoalSetAction>(
    node, raystarEndpoint(main_namespace_, "plan_goal_set"));
  ASSERT_TRUE(client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto goal = makeTestGoalSetActionGoal();
  // Both goals are much farther than this positive inclusive bound. The
  // complete result therefore contains zero paths without being a failed
  // planning operation.
  goal.max_path_lengths = {0.1, 0.1};
  goal.allow_self_crossing = true;
  goal.allow_unknown = true;
  auto goal_future = client->async_send_goal(goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto wrapped = result_future.get();

  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(wrapped.result, nullptr);
  EXPECT_TRUE(wrapped.result->success) << wrapped.result->message;
  EXPECT_EQ(wrapped.result->requested_start, goal.start);
  EXPECT_TRUE(wrapped.result->requested_allow_self_crossing);
  EXPECT_TRUE(wrapped.result->requested_allow_unknown);
  const auto& aggregate = wrapped.result->result_info;
  EXPECT_EQ(aggregate.status, PlanningResultInfo::STATUS_NO_PATH);
  EXPECT_TRUE(aggregate.request_satisfied);
  EXPECT_TRUE(aggregate.search_complete);
  EXPECT_TRUE(aggregate.output_complete);
  EXPECT_FALSE(aggregate.debug_requested);
  EXPECT_TRUE(aggregate.debug_output_complete);
  EXPECT_EQ(aggregate.requested_goal_count, 2u);
  EXPECT_EQ(aggregate.returned_goal_count, 2u);
  EXPECT_EQ(aggregate.completed_goal_count, 2u);
  EXPECT_EQ(aggregate.goals_with_paths, 0u);
  EXPECT_EQ(aggregate.found_path_count, 0u);
  EXPECT_EQ(aggregate.returned_path_count, 0u);
  ASSERT_EQ(wrapped.result->goal_results.size(), 2u);
  for (const auto& goal_result : wrapped.result->goal_results) {
    EXPECT_FALSE(goal_result.success);
    EXPECT_TRUE(goal_result.path_results.empty());
    EXPECT_EQ(goal_result.result_info.status, PlanningResultInfo::STATUS_NO_PATH);
    EXPECT_TRUE(goal_result.result_info.cost_bound_exhausted);
    EXPECT_TRUE(goal_result.result_info.search_complete);
    EXPECT_TRUE(goal_result.result_info.output_complete);
    EXPECT_TRUE(goal_result.result_info.request_satisfied);
    EXPECT_FALSE(goal_result.result_info.debug_requested);
    EXPECT_TRUE(goal_result.result_info.debug_output_complete);
  }
}

TEST_F(IntegrationTestFixture, GoalSetAllUnreachableStillUsesSucceededTransport) {
  const auto spawned = spawnRaystarNode("all_unreachable_goal_set");
  ASSERT_GT(spawned.pid, 0) << "Unable to exec all-unreachable raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  ChildProcessGuard process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_unreachable_goal_set_action_client");
  auto client = rclcpp_action::create_client<RaystarGoalSetAction>(
    node, raystarEndpoint(spawned.node_namespace, "plan_goal_set"));
  ASSERT_TRUE(client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto disconnected_map = makeTestGrid();
  std::fill(disconnected_map.data.begin(), disconnected_map.data.end(), 0);
  // Enclose both goal cells in a free pocket that is disconnected from start.
  for (unsigned int x = 18; x <= 26; ++x) {
    disconnected_map.data[8 * disconnected_map.info.width + x] = 100;
    disconnected_map.data[21 * disconnected_map.info.width + x] = 100;
  }
  for (unsigned int y = 8; y <= 21; ++y) {
    disconnected_map.data[y * disconnected_map.info.width + 18] = 100;
    disconnected_map.data[y * disconnected_map.info.width + 26] = 100;
  }
  ASSERT_TRUE(cacheMapAndWait(executor, node, disconnected_map, spawned.node_namespace, 5s));

  auto goal = makeTestGoalSetActionGoal();
  goal.map_id = raystar_interfaces::computeMapId(disconnected_map);
  goal.max_path_lengths = {100.0, 100.0};
  auto goal_future = client->async_send_goal(goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto wrapped = result_future.get();

  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(wrapped.result, nullptr);
  EXPECT_TRUE(wrapped.result->success) << wrapped.result->message;
  EXPECT_EQ(wrapped.result->result_info.status, PlanningResultInfo::STATUS_NO_PATH);
  EXPECT_TRUE(wrapped.result->result_info.request_satisfied);
  EXPECT_TRUE(wrapped.result->result_info.search_complete);
  EXPECT_TRUE(wrapped.result->result_info.output_complete);
  EXPECT_EQ(wrapped.result->result_info.completed_goal_count, 2u);
  EXPECT_EQ(wrapped.result->result_info.goals_with_paths, 0u);
  ASSERT_EQ(wrapped.result->goal_results.size(), 2u);
  for (const auto& goal_result : wrapped.result->goal_results) {
    EXPECT_FALSE(goal_result.success);
    EXPECT_TRUE(goal_result.path_results.empty());
    EXPECT_EQ(goal_result.result_info.status, PlanningResultInfo::STATUS_NO_PATH);
    EXPECT_TRUE(goal_result.result_info.cost_bound_exhausted);
    EXPECT_TRUE(goal_result.result_info.request_satisfied);
  }
}

TEST_F(IntegrationTestFixture, TransitionGraphActionReturnsCertifiedInteriorCrossingUpsPaths) {
  auto node = std::make_shared<rclcpp::Node>("test_transition_graph_action_client");
  auto planning_client = rclcpp_action::create_client<RaystarAction>(
    node, raystarEndpoint(main_namespace_, "plan_paths"));
  auto transition_client = rclcpp_action::create_client<RaystarTransitionAction>(
    node, raystarEndpoint(main_namespace_, "build_transition_graph"));
  ASSERT_TRUE(planning_client->wait_for_action_server(5s));
  ASSERT_TRUE(transition_client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto planning_goal = makeTestActionGoal();
  planning_goal.k = 2;
  auto planning_goal_future = planning_client->async_send_goal(planning_goal);
  ASSERT_TRUE(waitForFuture(executor, planning_goal_future, 5s));
  const auto planning_handle = planning_goal_future.get();
  ASSERT_NE(planning_handle, nullptr);
  auto planning_result_future = planning_client->async_get_result(planning_handle);
  ASSERT_TRUE(waitForFuture(executor, planning_result_future, 10s));
  const auto planning_wrapped = planning_result_future.get();
  ASSERT_EQ(planning_wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(planning_wrapped.result, nullptr);
  ASSERT_EQ(planning_wrapped.result->path_results.size(), 2u) << planning_wrapped.result->message;

  RaystarTransitionAction::Goal transition_goal;
  transition_goal.map_id = planning_goal.map_id;
  transition_goal.expected_environment_id = planning_wrapped.result->result_info.environment_id;
  transition_goal.allow_unknown = false;
  for (const auto& path_result : planning_wrapped.result->path_results) {
    ASSERT_GE(path_result.topology_path.poses.size(), 2u);
    transition_goal.tether_configurations.push_back(path_result.topology_path);
  }
  raystar_interfaces::msg::ConfigurationTransitionPair identity;
  identity.from_configuration = 0;
  identity.to_configuration = 0;
  raystar_interfaces::msg::ConfigurationTransitionPair change;
  change.from_configuration = 0;
  change.to_configuration = 1;
  transition_goal.transition_pairs = {identity, change};

  std::vector<RaystarTransitionAction::Feedback::ConstSharedPtr> feedback_messages;
  rclcpp_action::Client<RaystarTransitionAction>::SendGoalOptions transition_options;
  transition_options.feedback_callback =
    [&](auto, RaystarTransitionAction::Feedback::ConstSharedPtr feedback) {
      if (feedback)
        feedback_messages.push_back(std::move(feedback));
    };
  auto transition_goal_future =
    transition_client->async_send_goal(transition_goal, transition_options);
  ASSERT_TRUE(waitForFuture(executor, transition_goal_future, 5s));
  const auto transition_handle = transition_goal_future.get();
  ASSERT_NE(transition_handle, nullptr);
  auto transition_result_future = transition_client->async_get_result(transition_handle);
  ASSERT_TRUE(waitForFuture(executor, transition_result_future, 10s));
  const auto wrapped = transition_result_future.get();

  const auto feedback_deadline = std::chrono::steady_clock::now() + 2s;
  while ((feedback_messages.empty() || feedback_messages.back()->completed_transition_count != 2u ||
          feedback_messages.back()->stage != "transition batch complete") &&
         std::chrono::steady_clock::now() < feedback_deadline) {
    executor.spin_some(20ms);
  }

  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(wrapped.result, nullptr);
  EXPECT_TRUE(wrapped.result->success) << wrapped.result->message;
  EXPECT_EQ(wrapped.result->status, RaystarTransitionAction::Result::STATUS_COMPLETE);
  EXPECT_TRUE(raystar_interfaces::environmentIdsEqual(wrapped.result->environment_id,
                                                      transition_goal.expected_environment_id));
  EXPECT_EQ(wrapped.result->requested_transition_count, 2u);
  EXPECT_EQ(wrapped.result->completed_transition_count, 2u);
  ASSERT_EQ(wrapped.result->transitions.size(), 2u);
  EXPECT_LE(feedback_messages.size(), 100u);
  // ROS Action feedback is a volatile, non-authoritative live stream, not an
  // event log. This describes the API contract, not the RELIABLE QoS policy.
  // In particular, rclcpp discards feedback received before the accepted goal
  // response has installed the UUID-to-client-handle association. Therefore
  // an integration client can observe any ordered suffix/subsequence of the
  // producer trace. The reporter unit test verifies the complete producer
  // sequence; here we verify the transport-visible subsequence and use Result
  // above as the authoritative terminal progress record.
  std::uint32_t previous_completed = 0;
  int previous_stage_rank = -1;
  const auto stage_rank = [](const std::string& stage) {
    if (stage == "validating transition request")
      return 0;
    if (stage == "preparing transition environment")
      return 1;
    if (stage == "validating tether configurations")
      return 2;
    if (stage == "shortening transition pairs")
      return 3;
    if (stage == "transition batch complete")
      return 4;
    return -1;
  };
  for (const auto& feedback : feedback_messages) {
    ASSERT_NE(feedback, nullptr);
    EXPECT_EQ(feedback->requested_transition_count, 2u);
    EXPECT_GE(feedback->completed_transition_count, previous_completed);
    EXPECT_LE(feedback->completed_transition_count, 2u);
    EXPECT_FALSE(feedback->stage.empty());
    previous_completed = feedback->completed_transition_count;
    const int observed_stage_rank = stage_rank(feedback->stage);
    EXPECT_GE(observed_stage_rank, 0) << feedback->stage;
    EXPECT_GE(observed_stage_rank, previous_stage_rank) << feedback->stage;
    previous_stage_rank = observed_stage_rank;
  }
  EXPECT_DOUBLE_EQ(wrapped.result->transitions[0].path_length, 0.0);
  EXPECT_EQ(wrapped.result->transitions[0].path.poses.size(), 1u);
  EXPECT_GT(wrapped.result->transitions[1].path_length, 0.0);
  const auto& winding_transition = wrapped.result->transitions[1];
  bool winding_repeats_triangle_occurrence = false;
  for (size_t first = 0; first < winding_transition.triangle_occurrences.size(); ++first) {
    for (size_t second = first + 1; second < winding_transition.triangle_occurrences.size();
         ++second) {
      winding_repeats_triangle_occurrence =
        winding_repeats_triangle_occurrence || winding_transition.triangle_occurrences[first] ==
                                                 winding_transition.triangle_occurrences[second];
    }
  }
  EXPECT_TRUE(winding_repeats_triangle_occurrence)
    << "The ROS UPS certificate must retain lifted winding occurrences";
  EXPECT_NE(winding_transition.message.find("independently retraced"), std::string::npos)
    << winding_transition.message;
  for (const auto& transition : wrapped.result->transitions) {
    SCOPED_TRACE(transition.message);
    EXPECT_EQ(transition.status, raystar_interfaces::msg::HomotopyTransitionResult::STATUS_SUCCESS);
    EXPECT_TRUE(transition.collision_free);
    EXPECT_TRUE(transition.homotopy_preserved);
    EXPECT_TRUE(transition.locally_shortest);
    EXPECT_FALSE(transition.triangle_occurrences.empty());
    if (transition.path.poses.size() >= 2) {
      expectIndependentCollisionFreePath(makeTestGrid(), transition.path, false);
      raystar::ConservativeBinary64PathLength length_certificate;
      for (std::size_t pose_index = 1; pose_index < transition.path.poses.size(); ++pose_index) {
        const auto& first = transition.path.poses[pose_index - 1].pose.position;
        const auto& second = transition.path.poses[pose_index].pose.position;
        ASSERT_TRUE(length_certificate.addSegment(first.x, first.y, second.x, second.y));
      }
      double certified_length = std::numeric_limits<double>::quiet_NaN();
      ASSERT_TRUE(length_certificate.upperBound(certified_length));
      EXPECT_DOUBLE_EQ(transition.path_length, certified_length);
    } else {
      ASSERT_EQ(transition.path.poses.size(), 1u);
      EXPECT_DOUBLE_EQ(transition.path_length, 0.0);
    }
  }

  auto mismatched_goal = transition_goal;
  mismatched_goal.expected_environment_id.uuid[0] ^= 0xffU;
  auto mismatched_goal_future = transition_client->async_send_goal(mismatched_goal);
  ASSERT_TRUE(waitForFuture(executor, mismatched_goal_future, 5s));
  const auto mismatched_handle = mismatched_goal_future.get();
  ASSERT_NE(mismatched_handle, nullptr);
  auto mismatched_result_future = transition_client->async_get_result(mismatched_handle);
  ASSERT_TRUE(waitForFuture(executor, mismatched_result_future, 5s));
  const auto mismatched_wrapped = mismatched_result_future.get();
  EXPECT_EQ(mismatched_wrapped.code, rclcpp_action::ResultCode::ABORTED);
  ASSERT_NE(mismatched_wrapped.result, nullptr);
  EXPECT_FALSE(mismatched_wrapped.result->success);
  EXPECT_EQ(mismatched_wrapped.result->status,
            RaystarTransitionAction::Result::STATUS_INVALID_REQUEST);
  EXPECT_TRUE(raystar_interfaces::environmentIdsEqual(mismatched_wrapped.result->environment_id,
                                                      transition_goal.expected_environment_id));
  EXPECT_NE(mismatched_wrapped.result->message.find("expected_environment_id"), std::string::npos);
}

TEST_F(IntegrationTestFixture, TransitionGraphRejectsObstacleCrossingSharedPrefix) {
  auto node = std::make_shared<rclcpp::Node>("test_transition_shared_prefix_client");
  auto client = rclcpp_action::create_client<RaystarTransitionAction>(
    node, raystarEndpoint(main_namespace_, "build_transition_graph"));
  ASSERT_TRUE(client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  const auto pose = [](double x, double y) {
    geometry_msgs::msg::PoseStamped result;
    result.header.frame_id = "map";
    result.pose.position.x = x;
    result.pose.position.y = y;
    result.pose.orientation.w = 1.0;
    return result;
  };
  const auto path = [&](double endpoint_y) {
    nav_msgs::msg::Path result;
    result.header.frame_id = "map";
    // The shared 5.5,15.5 -> 25.5,15.5 segment crosses the central
    // obstacle. Both endpoints and the post-prefix pairwise reference are in
    // free space, so validating only alpha_a^{-1} * alpha_b would miss it.
    result.poses = {pose(5.5, 15.5), pose(25.5, 15.5), pose(25.5, endpoint_y)};
    return result;
  };

  RaystarTransitionAction::Goal goal;
  goal.map_id = raystar_interfaces::computeMapId(makeTestGrid());
  goal.allow_unknown = false;
  goal.tether_configurations = {path(14.5), path(16.5)};
  raystar_interfaces::msg::ConfigurationTransitionPair pair;
  pair.from_configuration = 0;
  pair.to_configuration = 1;
  goal.transition_pairs = {pair};

  auto goal_future = client->async_send_goal(goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto wrapped = result_future.get();

  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::ABORTED);
  ASSERT_NE(wrapped.result, nullptr);
  EXPECT_FALSE(wrapped.result->success);
  EXPECT_EQ(wrapped.result->status, RaystarTransitionAction::Result::STATUS_INVALID_REQUEST);
  EXPECT_EQ(wrapped.result->completed_transition_count, 0u);
  EXPECT_TRUE(wrapped.result->transitions.empty());
  EXPECT_NE(wrapped.result->message.find("Tether configuration 0"), std::string::npos);
  EXPECT_NE(wrapped.result->message.find("collision-free reference"), std::string::npos);
}

TEST_F(IntegrationTestFixture, TransitionGraphAcceptsDegenerateHomeIdentityReference) {
  auto node = std::make_shared<rclcpp::Node>("test_transition_home_identity_client");
  auto client = rclcpp_action::create_client<RaystarTransitionAction>(
    node, raystarEndpoint(main_namespace_, "build_transition_graph"));
  ASSERT_TRUE(client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  geometry_msgs::msg::PoseStamped base;
  base.header.frame_id = "map";
  base.pose.position.x = 5.5;
  base.pose.position.y = 15.5;
  base.pose.orientation.w = 1.0;
  nav_msgs::msg::Path home;
  home.header.frame_id = "map";
  home.poses = {base};

  RaystarTransitionAction::Goal goal;
  goal.map_id = raystar_interfaces::computeMapId(makeTestGrid());
  goal.allow_unknown = false;
  goal.tether_configurations = {home};
  raystar_interfaces::msg::ConfigurationTransitionPair identity;
  identity.from_configuration = 0;
  identity.to_configuration = 0;
  goal.transition_pairs = {identity};

  auto goal_future = client->async_send_goal(goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto wrapped = result_future.get();

  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(wrapped.result, nullptr);
  EXPECT_TRUE(wrapped.result->success) << wrapped.result->message;
  EXPECT_EQ(wrapped.result->status, RaystarTransitionAction::Result::STATUS_COMPLETE);
  ASSERT_EQ(wrapped.result->transitions.size(), 1u);
  const auto& transition = wrapped.result->transitions[0];
  EXPECT_EQ(transition.status, raystar_interfaces::msg::HomotopyTransitionResult::STATUS_SUCCESS);
  EXPECT_DOUBLE_EQ(transition.path_length, 0.0);
  ASSERT_EQ(transition.path.poses.size(), 1u);
  EXPECT_DOUBLE_EQ(transition.path.poses[0].pose.position.x, 5.5);
  EXPECT_DOUBLE_EQ(transition.path.poses[0].pose.position.y, 15.5);
  EXPECT_TRUE(transition.collision_free);
  EXPECT_TRUE(transition.homotopy_preserved);
  EXPECT_TRUE(transition.locally_shortest);
}

TEST_F(IntegrationTestFixture, TransitionPathLengthCoversSerializedAffineRounding) {
  const auto spawned = spawnRaystarNode("transition_affine_rounding");
  ASSERT_GT(spawned.pid, 0) << "Unable to exec affine-rounding raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  ChildProcessGuard process(spawned.pid);
  auto node = std::make_shared<rclcpp::Node>("test_transition_affine_rounding_client");
  auto client = rclcpp_action::create_client<RaystarTransitionAction>(
    node, raystarEndpoint(spawned.node_namespace, "build_transition_graph"));
  ASSERT_TRUE(client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto map = makeTestGrid();
  map.info.width = 12;
  map.info.height = 12;
  map.info.resolution = 0x1.99999ap-5f;
  map.info.origin.position.x = 0.3;
  map.info.origin.position.y = 0.0;
  map.data.assign(12U * 12U, 0);
  ASSERT_TRUE(cacheMapAndWait(executor, node, map, spawned.node_namespace));

  const auto pose_at = [&](double grid_x, double grid_y) {
    geometry_msgs::msg::PoseStamped pose;
    pose.header.frame_id = "map";
    pose.pose.position.x =
      std::fma(grid_x, static_cast<double>(map.info.resolution), map.info.origin.position.x);
    pose.pose.position.y =
      std::fma(grid_y, static_cast<double>(map.info.resolution), map.info.origin.position.y);
    pose.pose.orientation.w = 1.0;
    return pose;
  };
  const auto base = pose_at(2.0, 5.5);
  const auto first_endpoint = pose_at(3.0, 5.5);
  const auto second_endpoint = pose_at(4.0, 5.5);
  nav_msgs::msg::Path first_configuration;
  first_configuration.header.frame_id = "map";
  first_configuration.poses = {base, first_endpoint};
  nav_msgs::msg::Path second_configuration;
  second_configuration.header.frame_id = "map";
  second_configuration.poses = {base, second_endpoint};

  RaystarTransitionAction::Goal goal;
  goal.map_id = raystar_interfaces::computeMapId(map);
  goal.tether_configurations = {first_configuration, second_configuration};
  raystar_interfaces::msg::ConfigurationTransitionPair pair;
  pair.from_configuration = 0;
  pair.to_configuration = 1;
  goal.transition_pairs = {pair};
  goal.allow_unknown = false;

  auto goal_future = client->async_send_goal(goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto wrapped = result_future.get();

  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(wrapped.result, nullptr);
  ASSERT_TRUE(wrapped.result->success) << wrapped.result->message;
  ASSERT_EQ(wrapped.result->transitions.size(), 1u);
  const auto& transition = wrapped.result->transitions.front();
  ASSERT_EQ(transition.status, raystar_interfaces::msg::HomotopyTransitionResult::STATUS_SUCCESS);
  ASSERT_EQ(transition.path.poses.size(), 2u);

  raystar::ConservativeBinary64PathLength certificate;
  const auto& first = transition.path.poses[0].pose.position;
  const auto& second = transition.path.poses[1].pose.position;
  ASSERT_TRUE(certificate.addSegment(first.x, first.y, second.x, second.y));
  double certified_length = std::numeric_limits<double>::quiet_NaN();
  ASSERT_TRUE(certificate.upperBound(certified_length));
  constexpr double expected_certified_length = 0x1.99999a0000008p-5;
  const double old_nominal_length = static_cast<double>(map.info.resolution);
  EXPECT_DOUBLE_EQ(certified_length, expected_certified_length);
  EXPECT_GT(certified_length, old_nominal_length);
  EXPECT_DOUBLE_EQ(transition.path_length, certified_length);
}

TEST_F(IntegrationTestFixture, TransitionGraphRejectsNonTautConfigurationReference) {
  auto node = std::make_shared<rclcpp::Node>("test_transition_non_taut_client");
  auto client = rclcpp_action::create_client<RaystarTransitionAction>(
    node, raystarEndpoint(main_namespace_, "build_transition_graph"));
  ASSERT_TRUE(client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  const auto pose = [](double x, double y) {
    geometry_msgs::msg::PoseStamped result;
    result.header.frame_id = "map";
    result.pose.position.x = x;
    result.pose.position.y = y;
    result.pose.orientation.w = 1.0;
    return result;
  };
  nav_msgs::msg::Path detour;
  detour.header.frame_id = "map";
  detour.poses = {pose(5.5, 5.5), pose(5.5, 8.5), pose(8.5, 8.5)};

  RaystarTransitionAction::Goal goal;
  goal.map_id = raystar_interfaces::computeMapId(makeTestGrid());
  goal.allow_unknown = false;
  goal.reference_path_policy =
    RaystarTransitionAction::Goal::REFERENCE_PATHS_MUST_BE_TAUT;
  goal.tether_configurations = {detour};
  raystar_interfaces::msg::ConfigurationTransitionPair identity;
  identity.from_configuration = 0;
  identity.to_configuration = 0;
  goal.transition_pairs = {identity};

  auto goal_future = client->async_send_goal(goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto wrapped = result_future.get();

  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::ABORTED);
  ASSERT_NE(wrapped.result, nullptr);
  EXPECT_FALSE(wrapped.result->success);
  EXPECT_EQ(wrapped.result->status, RaystarTransitionAction::Result::STATUS_INVALID_REQUEST);
  EXPECT_TRUE(wrapped.result->transitions.empty());
  EXPECT_NE(wrapped.result->message.find("locally shortest (taut)"), std::string::npos);
}

TEST_F(IntegrationTestFixture, TransitionGraphRejectsUnknownReferencePathPolicy) {
  auto node = std::make_shared<rclcpp::Node>("test_transition_unknown_reference_policy_client");
  auto client = rclcpp_action::create_client<RaystarTransitionAction>(
    node, raystarEndpoint(main_namespace_, "build_transition_graph"));
  ASSERT_TRUE(client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  geometry_msgs::msg::PoseStamped base;
  base.header.frame_id = "map";
  base.pose.position.x = 5.5;
  base.pose.position.y = 5.5;
  base.pose.orientation.w = 1.0;
  nav_msgs::msg::Path home;
  home.header.frame_id = "map";
  home.poses = {base};

  RaystarTransitionAction::Goal goal;
  goal.map_id = raystar_interfaces::computeMapId(makeTestGrid());
  goal.allow_unknown = false;
  goal.reference_path_policy = 255;
  goal.tether_configurations = {home};
  raystar_interfaces::msg::ConfigurationTransitionPair identity;
  identity.from_configuration = 0;
  identity.to_configuration = 0;
  goal.transition_pairs = {identity};

  auto goal_future = client->async_send_goal(goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto wrapped = result_future.get();

  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::ABORTED);
  ASSERT_NE(wrapped.result, nullptr);
  EXPECT_FALSE(wrapped.result->success);
  EXPECT_EQ(wrapped.result->status, RaystarTransitionAction::Result::STATUS_INVALID_REQUEST);
  EXPECT_TRUE(wrapped.result->transitions.empty());
  EXPECT_NE(wrapped.result->message.find("reference_path_policy"), std::string::npos);
}

TEST_F(IntegrationTestFixture, TransitionGraphMayShortenUntautConfigurationReferences) {
  auto node = std::make_shared<rclcpp::Node>("test_transition_untaut_reference_client");
  auto client = rclcpp_action::create_client<RaystarTransitionAction>(
    node, raystarEndpoint(main_namespace_, "build_transition_graph"));
  ASSERT_TRUE(client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  const auto pose = [](double x, double y) {
    geometry_msgs::msg::PoseStamped result;
    result.header.frame_id = "map";
    result.pose.position.x = x;
    result.pose.position.y = y;
    result.pose.orientation.w = 1.0;
    return result;
  };
  nav_msgs::msg::Path home;
  home.header.frame_id = "map";
  home.poses = {pose(5.5, 5.5)};
  nav_msgs::msg::Path detour;
  detour.header.frame_id = "map";
  detour.poses = {pose(5.5, 5.5), pose(5.5, 8.5), pose(8.5, 8.5)};

  RaystarTransitionAction::Goal goal;
  goal.map_id = raystar_interfaces::computeMapId(makeTestGrid());
  goal.allow_unknown = false;
  goal.reference_path_policy =
    RaystarTransitionAction::Goal::REFERENCE_PATHS_MAY_BE_UNTAUT;
  goal.tether_configurations = {home, detour};
  raystar_interfaces::msg::ConfigurationTransitionPair pair;
  pair.from_configuration = 0;
  pair.to_configuration = 1;
  goal.transition_pairs = {pair};

  auto goal_future = client->async_send_goal(goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto wrapped = result_future.get();

  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(wrapped.result, nullptr);
  ASSERT_TRUE(wrapped.result->success) << wrapped.result->message;
  EXPECT_EQ(wrapped.result->status, RaystarTransitionAction::Result::STATUS_COMPLETE);
  ASSERT_EQ(wrapped.result->transitions.size(), 1u);
  const auto& transition = wrapped.result->transitions.front();
  EXPECT_EQ(transition.status, raystar_interfaces::msg::HomotopyTransitionResult::STATUS_SUCCESS);
  ASSERT_EQ(transition.path.poses.size(), 2u);
  EXPECT_DOUBLE_EQ(transition.path.poses.front().pose.position.x, 5.5);
  EXPECT_DOUBLE_EQ(transition.path.poses.front().pose.position.y, 5.5);
  EXPECT_DOUBLE_EQ(transition.path.poses.back().pose.position.x, 8.5);
  EXPECT_DOUBLE_EQ(transition.path.poses.back().pose.position.y, 8.5);
  EXPECT_NEAR(transition.path_length, std::hypot(3.0, 3.0), 1.0e-12);
  EXPECT_TRUE(transition.collision_free);
  EXPECT_TRUE(transition.homotopy_preserved);
  EXPECT_TRUE(transition.locally_shortest);
}

TEST_F(IntegrationTestFixture, TransitionGraphUsesRawContourAfterSimplifiedPlannerRequest) {
  const auto spawned = spawnRaystarNode("raw_contour_after_planner");
  ASSERT_GT(spawned.pid, 0) << "Unable to exec raw-contour raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  ChildProcessGuard process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_raw_contour_after_planner_client");
  auto planning_client = rclcpp_action::create_client<RaystarAction>(
    node, raystarEndpoint(spawned.node_namespace, "plan_paths"));
  auto transition_client = rclcpp_action::create_client<RaystarTransitionAction>(
    node, raystarEndpoint(spawned.node_namespace, "build_transition_graph"));
  ASSERT_TRUE(planning_client->wait_for_action_server(5s));
  ASSERT_TRUE(transition_client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  const auto map = makeRawContourWedgeMap();
  ASSERT_TRUE(cacheMapAndWait(executor, node, map, spawned.node_namespace));
  const auto pose = [](double x, double y) {
    geometry_msgs::msg::PoseStamped result;
    result.header.frame_id = "map";
    result.pose.position.x = x;
    result.pose.position.y = y;
    result.pose.orientation.w = 1.0;
    return result;
  };
  const auto base = pose(2.5, 47.5);
  const auto endpoint = pose(46.5, 3.5);

  // Populate the ordinary planner's simplified-contour state first.  On this
  // exact geometry, test_polymap.cpp proves that simplification removes a
  // raw-free wedge containing the untaut reference constructed below.
  RaystarAction::Goal planning_goal;
  planning_goal.map_id = raystar_interfaces::computeMapId(map);
  planning_goal.start = base;
  planning_goal.goal = endpoint;
  planning_goal.search_mode = RaystarAction::Goal::SEARCH_MODE_TOP_K;
  planning_goal.k = 1;
  planning_goal.allow_self_crossing = false;
  planning_goal.allow_unknown = false;
  planning_goal.include_debug = false;
  auto planning_goal_future = planning_client->async_send_goal(planning_goal);
  ASSERT_TRUE(waitForFuture(executor, planning_goal_future, 5s));
  const auto planning_handle = planning_goal_future.get();
  ASSERT_NE(planning_handle, nullptr);
  auto planning_result_future = planning_client->async_get_result(planning_handle);
  ASSERT_TRUE(waitForFuture(executor, planning_result_future, 10s));
  const auto planning_wrapped = planning_result_future.get();
  ASSERT_EQ(planning_wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(planning_wrapped.result, nullptr);
  ASSERT_TRUE(planning_wrapped.result->success) << planning_wrapped.result->message;
  ASSERT_EQ(planning_wrapped.result->path_results.size(), 1u);
  expectIndependentCollisionFreePath(
    map, planning_wrapped.result->path_results.front().topology_path, false);

  nav_msgs::msg::Path home;
  home.header.frame_id = "map";
  home.poses = {base};
  nav_msgs::msg::Path raw_free_untaut_reference;
  raw_free_untaut_reference.header.frame_id = "map";
  raw_free_untaut_reference.poses = {
    base, pose(32.0, 44.0), pose(39.0, 43.0), pose(45.0, 27.0), pose(46.0, 21.0), endpoint};
  expectIndependentCollisionFreePath(map, raw_free_untaut_reference, false);

  RaystarTransitionAction::Goal transition_goal;
  transition_goal.map_id = planning_goal.map_id;
  transition_goal.expected_environment_id = planning_wrapped.result->result_info.environment_id;
  transition_goal.reference_path_policy =
    RaystarTransitionAction::Goal::REFERENCE_PATHS_MAY_BE_UNTAUT;
  transition_goal.tether_configurations = {home, raw_free_untaut_reference};
  raystar_interfaces::msg::ConfigurationTransitionPair pair;
  pair.from_configuration = 0;
  pair.to_configuration = 1;
  transition_goal.transition_pairs = {pair};
  transition_goal.allow_unknown = false;

  auto transition_goal_future = transition_client->async_send_goal(transition_goal);
  ASSERT_TRUE(waitForFuture(executor, transition_goal_future, 5s));
  const auto transition_handle = transition_goal_future.get();
  ASSERT_NE(transition_handle, nullptr);
  auto transition_result_future = transition_client->async_get_result(transition_handle);
  ASSERT_TRUE(waitForFuture(executor, transition_result_future, 10s));
  const auto first_wrapped = transition_result_future.get();
  ASSERT_EQ(first_wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(first_wrapped.result, nullptr);
  ASSERT_TRUE(first_wrapped.result->success) << first_wrapped.result->message;
  EXPECT_EQ(first_wrapped.result->status, RaystarTransitionAction::Result::STATUS_COMPLETE);
  EXPECT_TRUE(raystar_interfaces::environmentIdsEqual(first_wrapped.result->environment_id,
                                                      transition_goal.expected_environment_id));
  ASSERT_EQ(first_wrapped.result->transitions.size(), 1u);
  const auto& first_transition = first_wrapped.result->transitions.front();
  EXPECT_EQ(first_transition.status,
            raystar_interfaces::msg::HomotopyTransitionResult::STATUS_SUCCESS);
  EXPECT_TRUE(first_transition.collision_free);
  EXPECT_TRUE(first_transition.homotopy_preserved);
  EXPECT_TRUE(first_transition.locally_shortest);
  ASSERT_EQ(first_transition.path.poses.size(), 2u);
  EXPECT_DOUBLE_EQ(first_transition.path.poses.front().pose.position.x, base.pose.position.x);
  EXPECT_DOUBLE_EQ(first_transition.path.poses.front().pose.position.y, base.pose.position.y);
  EXPECT_DOUBLE_EQ(first_transition.path.poses.back().pose.position.x, endpoint.pose.position.x);
  EXPECT_DOUBLE_EQ(first_transition.path.poses.back().pose.position.y, endpoint.pose.position.y);
  expectIndependentCollisionFreePath(map, first_transition.path, false);
  raystar::ConservativeBinary64PathLength first_length_certificate;
  ASSERT_TRUE(first_length_certificate.addSegment(base.pose.position.x,
                                                  base.pose.position.y,
                                                  endpoint.pose.position.x,
                                                  endpoint.pose.position.y));
  double expected_path_length = std::numeric_limits<double>::quiet_NaN();
  ASSERT_TRUE(first_length_certificate.upperBound(expected_path_length));
  EXPECT_DOUBLE_EQ(first_transition.path_length, expected_path_length);

  // A second identical public request exercises the raw-environment cache-hit
  // path without asserting on logs, timings, or any private cache structure.
  auto repeat_goal_future = transition_client->async_send_goal(transition_goal);
  ASSERT_TRUE(waitForFuture(executor, repeat_goal_future, 5s));
  const auto repeat_handle = repeat_goal_future.get();
  ASSERT_NE(repeat_handle, nullptr);
  auto repeat_result_future = transition_client->async_get_result(repeat_handle);
  ASSERT_TRUE(waitForFuture(executor, repeat_result_future, 10s));
  const auto repeat_wrapped = repeat_result_future.get();
  ASSERT_EQ(repeat_wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(repeat_wrapped.result, nullptr);
  ASSERT_TRUE(repeat_wrapped.result->success) << repeat_wrapped.result->message;
  ASSERT_EQ(repeat_wrapped.result->transitions.size(), 1u);
  const auto& repeated_transition = repeat_wrapped.result->transitions.front();
  EXPECT_EQ(repeated_transition.status, first_transition.status);
  EXPECT_EQ(repeated_transition.collision_free, first_transition.collision_free);
  EXPECT_EQ(repeated_transition.homotopy_preserved, first_transition.homotopy_preserved);
  EXPECT_EQ(repeated_transition.locally_shortest, first_transition.locally_shortest);
  EXPECT_EQ(repeated_transition.triangle_occurrences, first_transition.triangle_occurrences);
  EXPECT_DOUBLE_EQ(repeated_transition.path_length, first_transition.path_length);
  ASSERT_EQ(repeated_transition.path.poses.size(), first_transition.path.poses.size());
  for (size_t index = 0; index < repeated_transition.path.poses.size(); ++index) {
    EXPECT_DOUBLE_EQ(repeated_transition.path.poses[index].pose.position.x,
                     first_transition.path.poses[index].pose.position.x);
    EXPECT_DOUBLE_EQ(repeated_transition.path.poses[index].pose.position.y,
                     first_transition.path.poses[index].pose.position.y);
  }

  // MAY_BE_UNTAUT relaxes tautness only.  It must still reject a reference
  // whose intermediate waypoint lies in the occupied right-hand border.
  auto colliding_goal = transition_goal;
  nav_msgs::msg::Path colliding_reference;
  colliding_reference.header.frame_id = "map";
  colliding_reference.poses = {base, pose(49.5, 25.5), endpoint};
  colliding_goal.tether_configurations = {home, colliding_reference};
  auto colliding_goal_future = transition_client->async_send_goal(colliding_goal);
  ASSERT_TRUE(waitForFuture(executor, colliding_goal_future, 5s));
  const auto colliding_handle = colliding_goal_future.get();
  ASSERT_NE(colliding_handle, nullptr);
  auto colliding_result_future = transition_client->async_get_result(colliding_handle);
  ASSERT_TRUE(waitForFuture(executor, colliding_result_future, 10s));
  const auto colliding_wrapped = colliding_result_future.get();
  EXPECT_EQ(colliding_wrapped.code, rclcpp_action::ResultCode::ABORTED);
  ASSERT_NE(colliding_wrapped.result, nullptr);
  EXPECT_FALSE(colliding_wrapped.result->success);
  EXPECT_EQ(colliding_wrapped.result->status,
            RaystarTransitionAction::Result::STATUS_INVALID_REQUEST);
  EXPECT_EQ(colliding_wrapped.result->completed_transition_count, 0u);
  EXPECT_TRUE(colliding_wrapped.result->transitions.empty());
  EXPECT_NE(colliding_wrapped.result->message.find("Tether configuration 1"), std::string::npos);
  EXPECT_NE(colliding_wrapped.result->message.find("collision-free reference"),
            std::string::npos);
}

TEST_F(IntegrationTestFixture, ExhaustedSearchReportsFewerPathsWithoutStringParsing) {
  auto node = std::make_shared<rclcpp::Node>("test_fewer_paths_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
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
  EXPECT_EQ(response->result_info.status, PlanningResultInfo::STATUS_FEWER_PATHS);
  EXPECT_EQ(response->result_info.requested_path_count, 3u);
  EXPECT_EQ(response->result_info.found_path_count, 1u);
  EXPECT_EQ(response->result_info.returned_path_count, 1u);
  EXPECT_TRUE(response->result_info.search_complete);
  EXPECT_TRUE(response->result_info.output_complete);
  EXPECT_FALSE(response->result_info.request_satisfied);
  expectStructuredResultInvariants(*response);
}

TEST_F(IntegrationTestFixture, CostBoundedEnumerationIsCertifiedAcrossServiceAndAction) {
  auto node = std::make_shared<rclcpp::Node>("test_cost_bounded_client");
  auto service_client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
  auto action_client = rclcpp_action::create_client<RaystarAction>(
    node, raystarEndpoint(main_namespace_, "plan_paths"));
  ASSERT_TRUE(service_client->wait_for_service(2s));
  ASSERT_TRUE(action_client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto baseline_request = makeTestRequest();
  baseline_request->k = 3;
  const auto baseline = callService(executor, service_client, baseline_request);
  ASSERT_NE(baseline, nullptr);
  ASSERT_EQ(baseline->path_results.size(), 2u) << baseline->message;
  const double inclusive_bound = baseline->path_results.back().cost;

  auto bounded_request = makeTestRequest();
  bounded_request->search_mode = RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH;
  bounded_request->k = 0;
  bounded_request->max_path_length = inclusive_bound;
  const auto bounded = callService(executor, service_client, bounded_request);
  ASSERT_NE(bounded, nullptr);
  ASSERT_TRUE(bounded->success) << bounded->message;
  ASSERT_EQ(bounded->path_results.size(), baseline->path_results.size());
  EXPECT_EQ(bounded->result_info.status, PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_EQ(bounded->result_info.search_mode,
            RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH);
  EXPECT_DOUBLE_EQ(bounded->result_info.requested_max_path_length, inclusive_bound);
  EXPECT_EQ(bounded->result_info.requested_path_count, 0u);
  EXPECT_TRUE(bounded->result_info.cost_bound_exhausted);
  EXPECT_TRUE(bounded->result_info.search_complete);
  EXPECT_TRUE(bounded->result_info.output_complete);
  EXPECT_TRUE(bounded->result_info.request_satisfied);
  for (size_t i = 0; i < bounded->path_results.size(); ++i) {
    EXPECT_DOUBLE_EQ(bounded->path_results[i].cost, baseline->path_results[i].cost);
    EXPECT_LE(bounded->path_results[i].cost, inclusive_bound);
  }
  expectStructuredResultInvariants(*bounded);

  auto action_goal = makeTestActionGoal();
  action_goal.search_mode = RaystarAction::Goal::SEARCH_MODE_ALL_WITHIN_LENGTH;
  action_goal.k = 0;
  action_goal.max_path_length = inclusive_bound;
  auto goal_future = action_client->async_send_goal(action_goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = action_client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto action_result = result_future.get();
  EXPECT_EQ(action_result.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(action_result.result, nullptr);
  EXPECT_TRUE(action_result.result->result_info.cost_bound_exhausted);
  EXPECT_TRUE(action_result.result->result_info.request_satisfied);
  ASSERT_EQ(action_result.result->path_results.size(), bounded->path_results.size());
  for (size_t i = 0; i < bounded->path_results.size(); ++i)
    EXPECT_DOUBLE_EQ(action_result.result->path_results[i].cost, bounded->path_results[i].cost);
  expectStructuredResultInvariants(*action_result.result);

  auto below_request = makeTestRequest();
  below_request->search_mode = RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH;
  below_request->k = 0;
  below_request->max_path_length = std::nextafter(baseline->path_results.front().cost, 0.0);
  const auto below = callService(executor, service_client, below_request);
  ASSERT_NE(below, nullptr);
  EXPECT_FALSE(below->success);
  EXPECT_EQ(below->result_info.status, PlanningResultInfo::STATUS_NO_PATH);
  EXPECT_TRUE(below->result_info.cost_bound_exhausted);
  EXPECT_TRUE(below->result_info.search_complete);
  EXPECT_TRUE(below->result_info.output_complete);
  EXPECT_TRUE(below->result_info.request_satisfied);
  EXPECT_TRUE(below->path_results.empty());
  expectStructuredResultInvariants(*below);
}

TEST_F(IntegrationTestFixture, RadicalSumBoundarySurvivesServiceAndSharedTreeAction) {
  const auto spawned = spawnRaystarNode("radical_sum_boundary", {"max_path_points:=6"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec tight-bound raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  ChildProcessGuard process(spawned.pid);
  auto node = std::make_shared<rclcpp::Node>("test_radical_sum_boundary_client");
  auto service_client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  auto goal_set_client = rclcpp_action::create_client<RaystarGoalSetAction>(
    node, raystarEndpoint(spawned.node_namespace, "plan_goal_set"));
  ASSERT_TRUE(service_client->wait_for_service(5s));
  ASSERT_TRUE(goal_set_client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto map = makeTestGrid();
  map.info.width = 60;
  map.info.height = 60;
  map.info.resolution = 1.0f;
  map.data.assign(60U * 60U, 0);
  map.data[30U * 60U + 20U] = 100;
  ASSERT_TRUE(cacheMapAndWait(executor, node, map, spawned.node_namespace));

  constexpr double inclusive_bound = 0x1.79fa384f9da53p+4;
  auto request = makeTestRequest();
  request->map = map;
  request->start.pose.position.x = 19.5;
  request->start.pose.position.y = 52.5;
  request->goal.pose.position.x = 20.5;
  request->goal.pose.position.y = 29.0;
  request->search_mode = RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH;
  request->k = 0;
  request->max_path_length = inclusive_bound;
  const auto service_response = callService(executor, service_client, request);
  ASSERT_NE(service_response, nullptr);
  ASSERT_TRUE(service_response->success) << service_response->message;
  EXPECT_EQ(service_response->result_info.status, PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_EQ(service_response->result_info.found_path_count, 1u);
  EXPECT_EQ(service_response->result_info.returned_path_count, 1u);
  EXPECT_TRUE(service_response->result_info.cost_bound_exhausted);
  EXPECT_TRUE(service_response->result_info.output_complete);
  EXPECT_TRUE(service_response->result_info.request_satisfied);
  ASSERT_EQ(service_response->path_results.size(), 1u);
  const auto& service_path = service_response->path_results.front();
  EXPECT_DOUBLE_EQ(service_path.cost, inclusive_bound);
  EXPECT_LE(service_path.cost, request->max_path_length);
  // The dense interpolation cannot fit this node's six-point aggregate
  // budget; bounded mode must deterministically publish topology geometry for
  // both fields instead of reporting partial output.
  ASSERT_EQ(service_path.path.poses.size(), 3u);
  EXPECT_EQ(service_path.path, service_path.topology_path);
  expectIndependentCollisionFreePath(map, service_path.path, false);
  expectStructuredResultInvariants(*service_response);

  RaystarGoalSetAction::Goal goal;
  goal.map_id = raystar_interfaces::computeMapId(map);
  goal.start = request->start;
  goal.goals = {request->goal};
  goal.max_path_lengths = {inclusive_bound};
  goal.allow_self_crossing = false;
  goal.allow_unknown = false;
  goal.include_debug = false;
  auto goal_future = goal_set_client->async_send_goal(goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = goal_set_client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto wrapped = result_future.get();
  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(wrapped.result, nullptr);
  EXPECT_TRUE(wrapped.result->success) << wrapped.result->message;
  EXPECT_TRUE(wrapped.result->result_info.request_satisfied);
  EXPECT_TRUE(wrapped.result->result_info.output_complete);
  EXPECT_EQ(wrapped.result->result_info.found_path_count, 1u);
  ASSERT_EQ(wrapped.result->goal_results.size(), 1u);
  const auto& goal_result = wrapped.result->goal_results.front();
  EXPECT_TRUE(goal_result.success) << goal_result.message;
  EXPECT_EQ(goal_result.result_info.found_path_count, 1u);
  EXPECT_EQ(goal_result.result_info.returned_path_count, 1u);
  EXPECT_TRUE(goal_result.result_info.request_satisfied);
  ASSERT_EQ(goal_result.path_results.size(), 1u);
  EXPECT_DOUBLE_EQ(goal_result.path_results.front().cost, inclusive_bound);
  EXPECT_EQ(goal_result.path_results.front().path, goal_result.path_results.front().topology_path);
}

TEST_F(IntegrationTestFixture, MetricSearchSupersetExtrasRemainCompleteAndResourceNeutral) {
  const auto spawned =
    spawnRaystarNode("metric_search_superset", {"max_path_points:=4", "max_cost_bounded_paths:=1"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec metric-superset raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  ChildProcessGuard process(spawned.pid);
  auto node = std::make_shared<rclcpp::Node>("test_metric_search_superset_client");
  auto service_client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  auto goal_set_client = rclcpp_action::create_client<RaystarGoalSetAction>(
    node, raystarEndpoint(spawned.node_namespace, "plan_goal_set"));
  ASSERT_TRUE(service_client->wait_for_service(5s));
  ASSERT_TRUE(goal_set_client->wait_for_action_server(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto map = makeTestGrid();
  map.info.width = 60;
  map.info.height = 60;
  map.info.resolution = 0x1.99999ap-5f;
  map.data.assign(60U * 60U, 0);
  ASSERT_TRUE(cacheMapAndWait(executor, node, map, spawned.node_namespace));
  constexpr double inclusive_bound = 0x1.cd82b4b9764c5p-4;
  const double excluded_bound = std::nextafter(inclusive_bound, 0.0);
  const double resolution = static_cast<double>(map.info.resolution);

  auto request = makeTestRequest();
  request->map = map;
  request->start.pose.position.x = 20.875 * resolution;
  request->start.pose.position.y = 10.125 * resolution;
  request->goal.pose.position.x = 20.75 * resolution;
  request->goal.pose.position.y = 12.375 * resolution;
  request->search_mode = RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH;
  request->k = 0;
  request->max_path_length = inclusive_bound;
  const auto accepted = callService(executor, service_client, request);
  ASSERT_NE(accepted, nullptr);
  ASSERT_TRUE(accepted->success) << accepted->message;
  ASSERT_EQ(accepted->path_results.size(), 1u);
  EXPECT_DOUBLE_EQ(accepted->path_results.front().cost, inclusive_bound);
  EXPECT_EQ(accepted->path_results.front().path, accepted->path_results.front().topology_path);
  EXPECT_TRUE(accepted->result_info.output_complete);
  EXPECT_TRUE(accepted->result_info.request_satisfied);

  request->max_path_length = excluded_bound;
  const auto excluded = callService(executor, service_client, request);
  ASSERT_NE(excluded, nullptr);
  EXPECT_FALSE(excluded->success);
  EXPECT_TRUE(excluded->path_results.empty());
  EXPECT_EQ(excluded->result_info.status, PlanningResultInfo::STATUS_NO_PATH);
  EXPECT_EQ(excluded->result_info.found_path_count, 0u);
  EXPECT_EQ(excluded->result_info.returned_path_count, 0u);
  EXPECT_EQ(excluded->result_info.limits_reached, PlanningResultInfo::LIMIT_NONE);
  EXPECT_TRUE(excluded->result_info.search_complete);
  EXPECT_TRUE(excluded->result_info.cost_bound_exhausted);
  EXPECT_TRUE(excluded->result_info.output_complete);
  EXPECT_TRUE(excluded->result_info.request_satisfied);
  expectStructuredResultInvariants(*excluded);

  RaystarGoalSetAction::Goal goal;
  goal.map_id = raystar_interfaces::computeMapId(map);
  goal.start = request->start;
  goal.goals = {request->goal};
  goal.max_path_lengths = {excluded_bound};
  goal.allow_self_crossing = false;
  goal.allow_unknown = false;
  goal.include_debug = false;
  auto goal_future = goal_set_client->async_send_goal(goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s));
  const auto goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr);
  auto result_future = goal_set_client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s));
  const auto wrapped = result_future.get();
  EXPECT_EQ(wrapped.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(wrapped.result, nullptr);
  EXPECT_TRUE(wrapped.result->success) << wrapped.result->message;
  EXPECT_EQ(wrapped.result->result_info.status, PlanningResultInfo::STATUS_NO_PATH);
  EXPECT_EQ(wrapped.result->result_info.found_path_count, 0u);
  EXPECT_EQ(wrapped.result->result_info.returned_path_count, 0u);
  EXPECT_TRUE(wrapped.result->result_info.search_complete);
  EXPECT_TRUE(wrapped.result->result_info.output_complete);
  EXPECT_TRUE(wrapped.result->result_info.request_satisfied);
  ASSERT_EQ(wrapped.result->goal_results.size(), 1u);
  const auto& empty_goal = wrapped.result->goal_results.front();
  EXPECT_FALSE(empty_goal.success);
  EXPECT_TRUE(empty_goal.path_results.empty());
  EXPECT_EQ(empty_goal.result_info.status, PlanningResultInfo::STATUS_NO_PATH);
  EXPECT_EQ(empty_goal.result_info.found_path_count, 0u);
  EXPECT_TRUE(empty_goal.result_info.output_complete);
  EXPECT_TRUE(empty_goal.result_info.request_satisfied);
}

TEST_F(IntegrationTestFixture, SerializedLengthCertificateFiltersRoundedDownSearchSuperset) {
  auto node = std::make_shared<rclcpp::Node>("test_exact_serialized_length_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto request = makeTestRequest();
  std::fill(request->map.data.begin(), request->map.data.end(), 0);
  request->start.pose.position.x = 5.25;
  request->start.pose.position.y = 5.25;
  request->goal.pose.position.x = 6.25;
  request->goal.pose.position.y = 5.25 + std::ldexp(1.0, -27);
  request->k = 1;

  // Exact squared length is 1 + 2^-54, although nearest-rounded hypot is 1.
  ASSERT_DOUBLE_EQ(std::hypot(request->goal.pose.position.x - request->start.pose.position.x,
                              request->goal.pose.position.y - request->start.pose.position.y),
                   1.0);
  const auto baseline = callService(executor, client, request);
  ASSERT_NE(baseline, nullptr);
  ASSERT_TRUE(baseline->success) << baseline->message;
  ASSERT_EQ(baseline->path_results.size(), 1u);
  EXPECT_DOUBLE_EQ(baseline->path_results.front().cost,
                   std::nextafter(1.0, std::numeric_limits<double>::infinity()));
  expectStructuredResultInvariants(*baseline);

  request->search_mode = RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH;
  request->k = 0;
  request->max_path_length = 1.0;
  const auto bounded = callService(executor, client, request);
  ASSERT_NE(bounded, nullptr);
  EXPECT_FALSE(bounded->success);
  EXPECT_TRUE(bounded->path_results.empty());
  // Padding deliberately keeps the Core search a superset, but the
  // authoritative serialized-world topology certificate classifies this
  // path outside the original inclusive metric request.  It is therefore
  // neither a found path nor an output omission.
  EXPECT_EQ(bounded->result_info.found_path_count, 0u);
  EXPECT_EQ(bounded->result_info.returned_path_count, 0u);
  EXPECT_EQ(bounded->result_info.status, PlanningResultInfo::STATUS_NO_PATH);
  EXPECT_TRUE(bounded->result_info.search_complete);
  EXPECT_TRUE(bounded->result_info.cost_bound_exhausted);
  EXPECT_TRUE(bounded->result_info.output_complete);
  EXPECT_TRUE(bounded->result_info.request_satisfied);
  expectStructuredResultInvariants(*bounded);
}

TEST_F(IntegrationTestFixture, RejectsAmbiguousCostBoundedRequestFieldsThenRecovers) {
  auto node = std::make_shared<rclcpp::Node>("test_cost_bounded_validation_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  const std::vector<std::function<void(RaystarService::Request&)>> invalid_mutations = {
    [](auto& request) {
      request.search_mode = RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH;
      request.k = 1;
      request.max_path_length = 100.0;
    },
    [](auto& request) {
      request.search_mode = RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH;
      request.k = 0;
      request.max_path_length = 0.0;
    },
    [](auto& request) { request.search_mode = 255; },
    [](auto& request) { request.max_path_length = 1.0; },
  };
  for (const auto& mutate : invalid_mutations) {
    auto request = makeTestRequest();
    mutate(*request);
    const auto response = callService(executor, client, request);
    ASSERT_NE(response, nullptr);
    EXPECT_FALSE(response->success);
    EXPECT_EQ(response->result_info.status, PlanningResultInfo::STATUS_INVALID_REQUEST);
    EXPECT_FALSE(response->result_info.cost_bound_exhausted);
    EXPECT_TRUE(response->path_results.empty());
  }

  const auto recovered = callService(executor, client, makeTestRequest());
  ASSERT_NE(recovered, nullptr);
  EXPECT_TRUE(recovered->success) << recovered->message;
}

TEST_F(IntegrationTestFixture, SelfCrossingPolicyIsEnforcedAcrossServiceAndActionOutputs) {
  auto node = std::make_shared<rclcpp::Node>("test_self_crossing_contract_client");
  auto service_client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
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
  const auto service_response = callService(executor, service_client, service_request);
  ASSERT_NE(service_response, nullptr);
  expectSingleObstacleSimplePathContract(service_request->map, *service_response);

  ASSERT_TRUE(cacheMapAndWait(executor, node, service_request->map, main_namespace_))
    << "Raystar server did not admit the self-crossing contract map";
  auto action_goal = makeTestActionGoal();
  action_goal.k = 3;
  action_goal.allow_self_crossing = false;
  auto goal_future = action_client->async_send_goal(action_goal);
  ASSERT_TRUE(waitForFuture(executor, goal_future, 5s))
    << "Self-crossing-disabled Action goal was not acknowledged";
  const RaystarGoalHandle::SharedPtr goal_handle = goal_future.get();
  ASSERT_NE(goal_handle, nullptr) << "Self-crossing-disabled Action goal was rejected";
  auto result_future = action_client->async_get_result(goal_handle);
  ASSERT_TRUE(waitForFuture(executor, result_future, 10s))
    << "Self-crossing-disabled Action goal did not complete";
  const auto action_result = result_future.get();
  EXPECT_EQ(action_result.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(action_result.result, nullptr);
  expectSingleObstacleSimplePathContract(service_request->map, *action_result.result);

  // Enabling the policy only permits self-crossing candidates; it does not
  // require the planner to manufacture one.  This smoke request therefore
  // verifies acceptance and the ordinary output contract without asserting a
  // particular self-intersecting path or a policy-dependent path count.
  auto permissive_goal = makeTestActionGoal();
  permissive_goal.k = 1;
  permissive_goal.allow_self_crossing = true;
  auto permissive_goal_future = action_client->async_send_goal(permissive_goal);
  ASSERT_TRUE(waitForFuture(executor, permissive_goal_future, 5s))
    << "Self-crossing-enabled Action goal was not acknowledged";
  const RaystarGoalHandle::SharedPtr permissive_handle = permissive_goal_future.get();
  ASSERT_NE(permissive_handle, nullptr) << "Self-crossing-enabled Action goal was rejected";
  auto permissive_result_future = action_client->async_get_result(permissive_handle);
  ASSERT_TRUE(waitForFuture(executor, permissive_result_future, 10s))
    << "Self-crossing-enabled Action goal did not complete";
  const auto permissive_result = permissive_result_future.get();
  EXPECT_EQ(permissive_result.code, rclcpp_action::ResultCode::SUCCEEDED);
  ASSERT_NE(permissive_result.result, nullptr);
  const auto& permissive_response = *permissive_result.result;
  EXPECT_TRUE(permissive_response.success) << permissive_response.message;
  EXPECT_EQ(permissive_response.result_info.status, PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_EQ(permissive_response.result_info.limits_reached, PlanningResultInfo::LIMIT_NONE);
  EXPECT_EQ(permissive_response.result_info.requested_path_count, 1u);
  EXPECT_EQ(permissive_response.result_info.found_path_count, 1u);
  EXPECT_EQ(permissive_response.result_info.returned_path_count, 1u);
  EXPECT_TRUE(permissive_response.result_info.request_satisfied);
  EXPECT_TRUE(permissive_response.result_info.search_complete);
  EXPECT_TRUE(permissive_response.result_info.output_complete);
  EXPECT_FALSE(permissive_response.result_info.debug_requested);
  EXPECT_TRUE(permissive_response.result_info.debug_output_complete);
  EXPECT_TRUE(
    raystar_interfaces::mapIdsEqual(permissive_response.result_info.map_id,
                                    raystar_interfaces::computeMapId(service_request->map)));
  ASSERT_EQ(permissive_response.path_results.size(), 1u);
  expectIndependentCollisionFreePath(
    service_request->map, permissive_response.path_results.front().path, false);
  expectStructuredResultInvariants(permissive_response);
}

TEST_F(IntegrationTestFixture, OccupancyThresholdAndUnknownPolicyAreApplied) {
  auto node = std::make_shared<rclcpp::Node>("test_occupancy_threshold_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto below_threshold = callService(executor, client, makeVerticalBarrierRequest(98));
  ASSERT_NE(below_threshold, nullptr);
  EXPECT_TRUE(below_threshold->success) << below_threshold->message;

  auto at_threshold = callService(executor, client, makeVerticalBarrierRequest(99));
  ASSERT_NE(at_threshold, nullptr);
  EXPECT_FALSE(at_threshold->success);
  EXPECT_NE(at_threshold->message.find("No path exists"), std::string::npos)
    << at_threshold->message;

  auto unknown_blocked = callService(executor, client, makeVerticalBarrierRequest(-1, false));
  ASSERT_NE(unknown_blocked, nullptr);
  EXPECT_FALSE(unknown_blocked->success);
  EXPECT_NE(unknown_blocked->message.find("No path exists"), std::string::npos)
    << unknown_blocked->message;

  auto unknown_allowed = callService(executor, client, makeVerticalBarrierRequest(-1, true));
  ASSERT_NE(unknown_allowed, nullptr);
  EXPECT_TRUE(unknown_allowed->success) << unknown_allowed->message;
}

TEST_F(IntegrationTestFixture, FixedSeedSmallGridDifferentialProperties) {
  auto node = std::make_shared<rclcpp::Node>("test_small_grid_property_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  constexpr size_t property_case_count = 24;
  for (size_t case_index = 0; case_index < property_case_count; ++case_index) {
    const auto property_case = raystar::test_property::makeSmallGridCase(case_index);
    SCOPED_TRACE(raystar::test_property::describeSmallGridCase(property_case));

    const bool bfs_reachable = raystar::test_property::fourNeighborReachable(property_case);
    ASSERT_EQ(bfs_reachable, property_case.expected_reachable)
      << "The fixed-seed generator did not preserve its reachability label";
    ASSERT_FALSE(raystar::test_property::hasUnsupportedDiagonalContact(property_case))
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
    if (bfs_reachable) {
      expectSmallGridReachableResult(property_case, *request, *first_response);
      expectSmallGridReachableResult(property_case, *request, *second_response);
    } else {
      expectSmallGridUnreachableResult(property_case, *first_response);
      expectSmallGridUnreachableResult(property_case, *second_response);
    }
    expectStableSmallGridResponses(*first_response, *second_response);
  }
}

TEST_F(IntegrationTestFixture, ActionCancelRejectsConcurrentGoalAndRecovers) {
  auto node = std::make_shared<rclcpp::Node>("test_action_cancel_client");
  auto client = rclcpp_action::create_client<RaystarAction>(
    node, raystarEndpoint(main_namespace_, "plan_paths"));
  auto service_client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_action_server(5s));
  ASSERT_TRUE(service_client->wait_for_service(5s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  ASSERT_TRUE(cacheMapAndWait(executor, node, makeLongRunningActionMap(), main_namespace_))
    << "The long-running map was not admitted by the Action server";
  auto bounded_long_goal = makeLongRunningActionGoal();
  bounded_long_goal.search_mode = RaystarAction::Goal::SEARCH_MODE_ALL_WITHIN_LENGTH;
  bounded_long_goal.k = 0;
  bounded_long_goal.max_path_length = 1.0e9;
  auto first_goal_future = client->async_send_goal(bounded_long_goal);
  ASSERT_TRUE(waitForFuture(executor, first_goal_future, 10s))
    << "The long-running Action goal was not acknowledged";
  const RaystarGoalHandle::SharedPtr first_handle = first_goal_future.get();
  ASSERT_NE(first_handle, nullptr) << "The first Action goal was rejected";
  auto first_result_future = client->async_get_result(first_handle);

  // Goal admission is serialized by the server.  Waiting for the first goal
  // handle before sending this request makes the assertion independent of
  // discovery/transport ordering: the first goal is already the active goal.
  auto second_goal_future = client->async_send_goal(makeTestActionGoal());
  const bool second_response_ready = waitForFuture(executor, second_goal_future, 5s);
  RaystarGoalHandle::SharedPtr second_handle;
  if (second_response_ready)
    second_handle = second_goal_future.get();
  else
    ADD_FAILURE() << "The concurrent Action goal received no admission response";

  auto busy_service_future = service_client->async_send_request(makeTestRequest());
  const bool busy_service_ready = waitForFuture(executor, busy_service_future, 5s);
  RaystarService::Response::SharedPtr busy_service_response;
  if (busy_service_ready)
    busy_service_response = busy_service_future.get();
  else
    ADD_FAILURE() << "The compatibility Service received no busy response";

  auto cancel_future = client->async_cancel_goal(first_handle);
  const bool cancel_response_ready = waitForFuture(executor, cancel_future, 5s);
  if (!cancel_response_ready) {
    ADD_FAILURE() << "Cancel request received no response";
  } else {
    const auto cancel_response = cancel_future.get();
    EXPECT_FALSE(cancel_response->goals_canceling.empty())
      << "The active goal was not accepted for cancellation";
  }

  // If a broken server accepted the second goal, request its cancellation as
  // cleanup before making the rejection assertion below.
  if (second_handle)
    (void)client->async_cancel_goal(second_handle);

  const bool canceled_result_ready = waitForFuture(executor, first_result_future, 15s);
  ASSERT_TRUE(canceled_result_ready) << "Canceled Action goal did not reach a terminal state";
  const auto canceled_result = first_result_future.get();

  EXPECT_TRUE(second_response_ready);
  EXPECT_EQ(second_handle, nullptr)
    << "A second Action goal was accepted while planning was active";
  EXPECT_TRUE(busy_service_ready);
  ASSERT_NE(busy_service_response, nullptr);
  EXPECT_FALSE(busy_service_response->success);
  EXPECT_EQ(busy_service_response->result_info.status, PlanningResultInfo::STATUS_BUSY);
  EXPECT_EQ(busy_service_response->result_info.requested_path_count, 3u);
  expectStructuredResultInvariants(*busy_service_response);
  EXPECT_NE(lowercase(busy_service_response->message).find("busy"), std::string::npos)
    << busy_service_response->message;
  EXPECT_EQ(canceled_result.code, rclcpp_action::ResultCode::CANCELED);
  ASSERT_NE(canceled_result.result, nullptr);
  EXPECT_FALSE(canceled_result.result->success);
  EXPECT_EQ(canceled_result.result->result_info.status, PlanningResultInfo::STATUS_CANCELLED);
  EXPECT_NE(
    canceled_result.result->result_info.limits_reached & PlanningResultInfo::LIMIT_CANCELLED, 0u);
  EXPECT_FALSE(canceled_result.result->result_info.search_complete);
  EXPECT_FALSE(canceled_result.result->result_info.output_complete);
  EXPECT_FALSE(canceled_result.result->result_info.request_satisfied);
  EXPECT_FALSE(canceled_result.result->result_info.cost_bound_exhausted);
  EXPECT_EQ(canceled_result.result->result_info.search_mode,
            RaystarAction::Goal::SEARCH_MODE_ALL_WITHIN_LENGTH);
  EXPECT_DOUBLE_EQ(canceled_result.result->result_info.requested_max_path_length,
                   bounded_long_goal.max_path_length);
  EXPECT_EQ(canceled_result.result->result_info.requested_path_count, 0u);
  EXPECT_TRUE(canceled_result.result->path_results.empty());
  expectStructuredResultInvariants(*canceled_result.result);
  EXPECT_NE(lowercase(canceled_result.result->message).find("cancel"), std::string::npos)
    << canceled_result.result->message;

  ASSERT_TRUE(cacheMapAndWait(executor, node, makeTestGrid(), main_namespace_))
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
  EXPECT_TRUE(recovery_result.result->success) << recovery_result.result->message;
  EXPECT_FALSE(recovery_result.result->path_results.empty());
  EXPECT_TRUE(recovery_result.result->result_info.search_complete);
  EXPECT_TRUE(recovery_result.result->result_info.output_complete);
  expectStructuredResultInvariants(*recovery_result.result);
}

TEST_F(IntegrationTestFixture, SigintDuringActiveActionJoinsWorkerAndExits) {
  const auto spawned = spawnRaystarNode("graceful_shutdown", {"planning_timeout_ms:=60000"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec graceful-shutdown raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  const pid_t shutdown_pid = spawned.pid;
  ChildProcessGuard shutdown_process(shutdown_pid);

  auto node = std::make_shared<rclcpp::Node>("test_action_shutdown_client");
  auto client = rclcpp_action::create_client<RaystarAction>(
    node, raystarEndpoint(spawned.node_namespace, "plan_paths"));
  ASSERT_TRUE(client->wait_for_action_server(5s));

  bool worker_started = false;
  std::size_t deleteall_count = 0;
  auto worker_subscription = node->create_subscription<visualization_msgs::msg::MarkerArray>(
    raystarEndpoint(spawned.node_namespace, "non_homotopic_paths"),
    rclcpp::QoS(rclcpp::KeepLast(10)).transient_local(),
    [&](const visualization_msgs::msg::MarkerArray::ConstSharedPtr message) {
      const auto count = static_cast<std::size_t>(
        std::count_if(message->markers.begin(), message->markers.end(), [](const auto& marker) {
          return marker.action == visualization_msgs::msg::Marker::DELETEALL;
        }));
      deleteall_count += count;
      worker_started = deleteall_count > 0;
    });
  (void)worker_subscription;

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);
  ASSERT_TRUE(cacheMapAndWait(executor, node, makeLongRunningActionMap(), spawned.node_namespace))
    << "The shutdown test map was not admitted";
  auto goal_future = client->async_send_goal(makeLongRunningActionGoal());
  ASSERT_TRUE(waitForFuture(executor, goal_future, 10s))
    << "The shutdown test Action goal was not acknowledged";
  ASSERT_NE(goal_future.get(), nullptr) << "The shutdown test Action goal was rejected";

  // The Action handle proves admission; the first replacement marker proves
  // the owned worker has entered executePlanning().  Shutdown is therefore
  // synchronized to an active worker rather than an arbitrary delay.
  const auto worker_deadline = std::chrono::steady_clock::now() + 5s;
  while (!worker_started && std::chrono::steady_clock::now() < worker_deadline) {
    executor.spin_some(20ms);
    std::this_thread::sleep_for(1ms);
  }
  ASSERT_TRUE(worker_started) << "The Action worker did not enter planning before shutdown";

  ASSERT_EQ(kill(shutdown_pid, SIGINT), 0);

  int child_status = 0;
  pid_t wait_result = 0;
  bool wait_failed = false;
  const auto exit_deadline = std::chrono::steady_clock::now() + 10s;
  while (std::chrono::steady_clock::now() < exit_deadline) {
    wait_result = waitpid(shutdown_pid, &child_status, WNOHANG);
    if (wait_result == shutdown_pid)
      break;
    if (wait_result < 0) {
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
  while (deleteall_count < 2 && std::chrono::steady_clock::now() < clear_deadline) {
    executor.spin_some(20ms);
    std::this_thread::sleep_for(1ms);
  }
  // The first DELETEALL is the replacement marker that proves the worker was
  // active.  A second snapshot is published by the child destructor, but DDS
  // may drop an in-flight transient-local sample as the process exits.  Keep
  // this cross-process assertion at the observable best-effort boundary.
  EXPECT_GE(deleteall_count, 1u) << "The active worker did not publish a DELETEALL snapshot";
}

TEST_F(IntegrationTestFixture, PreservesContinuousEndpointsAndMetricCost) {
  auto node = std::make_shared<rclcpp::Node>("test_continuous_endpoint_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto request = makeTestRequest();
  request->map.info.resolution = 0.1f;
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
  for (size_t i = 1; i < path.poses.size(); ++i) {
    const auto& previous = path.poses[i - 1].pose.position;
    const auto& current = path.poses[i].pose.position;
    sampled_length += std::hypot(current.x - previous.x, current.y - previous.y);
  }
  EXPECT_NEAR(path_result.cost, sampled_length, 1e-8);
  EXPECT_GE(
    path_result.cost,
    std::hypot(goal_world.first - start_world.first, goal_world.second - start_world.second));
  expectIndependentCollisionFreePath(request->map, path, request->allow_unknown);
  EXPECT_EQ(response->result_info.status, PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_TRUE(response->result_info.request_satisfied);
  expectStructuredResultInvariants(*response);

  // Reusing the public metric cost as an inclusive bound must retain the
  // same shortest path even when metre/grid conversion uses a non-unit,
  // non-binary-exact float resolution.
  request->search_mode = RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH;
  request->k = 0;
  request->max_path_length = path_result.cost;
  const auto bounded = callService(executor, client, request);
  ASSERT_NE(bounded, nullptr);
  ASSERT_TRUE(bounded->success) << bounded->message;
  ASSERT_FALSE(bounded->path_results.empty());
  EXPECT_DOUBLE_EQ(bounded->path_results.front().cost, path_result.cost);
  for (const auto& bounded_path : bounded->path_results) {
    EXPECT_LE(bounded_path.cost, request->max_path_length);
    for (const auto* serialized_path : {&bounded_path.path, &bounded_path.topology_path}) {
      long double serialized_length = 0.0L;
      for (size_t index = 1; index < serialized_path->poses.size(); ++index) {
        const auto& left = serialized_path->poses[index - 1].pose.position;
        const auto& right = serialized_path->poses[index].pose.position;
        serialized_length +=
          static_cast<long double>(std::hypot(right.x - left.x, right.y - left.y));
      }
      EXPECT_LE(serialized_length, static_cast<long double>(request->max_path_length));
    }
  }
  EXPECT_TRUE(bounded->result_info.cost_bound_exhausted);
  EXPECT_TRUE(bounded->result_info.output_complete);
  EXPECT_TRUE(bounded->result_info.request_satisfied);
  expectStructuredResultInvariants(*bounded);

  request->max_path_length = std::nextafter(path_result.cost, 0.0);
  const auto just_below = callService(executor, client, request);
  ASSERT_NE(just_below, nullptr);
  EXPECT_FALSE(just_below->success);
  EXPECT_TRUE(just_below->path_results.empty());
  EXPECT_EQ(just_below->result_info.status, PlanningResultInfo::STATUS_NO_PATH);
  EXPECT_EQ(just_below->result_info.found_path_count, 0u);
  EXPECT_EQ(just_below->result_info.returned_path_count, 0u);
  EXPECT_TRUE(just_below->result_info.search_complete);
  EXPECT_TRUE(just_below->result_info.cost_bound_exhausted);
  EXPECT_TRUE(just_below->result_info.output_complete);
  EXPECT_TRUE(just_below->result_info.request_satisfied);
  expectStructuredResultInvariants(*just_below);
}

TEST_F(IntegrationTestFixture, RejectsOccupiedAndForcedBorderEndpointsSymmetricallyThenRecovers) {
  auto node = std::make_shared<rclcpp::Node>("test_endpoint_policy_matrix_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  using Mutator = std::function<void(RaystarService::Request&)>;
  struct InvalidEndpointCase {
    const char* name;
    const char* expected_endpoint;
    const char* expected_reason;
    Mutator mutate;
  };
  const std::vector<InvalidEndpointCase> invalid_cases = {
    {"occupied start",
     "invalid start",
     "occupied",
     [](auto& request) {
       request.start.pose.position.x = 15.5;
       request.start.pose.position.y = 15.5;
     }},
    {"occupied goal",
     "invalid goal",
     "occupied",
     [](auto& request) {
       request.goal.pose.position.x = 15.5;
       request.goal.pose.position.y = 15.5;
     }},
    {"forced-border start",
     "invalid start",
     "map boundary",
     [](auto& request) {
       request.start.pose.position.x = 0.5;
       request.start.pose.position.y = 15.5;
     }},
    {"forced-border goal",
     "invalid goal",
     "map boundary",
     [](auto& request) {
       request.goal.pose.position.x = 29.5;
       request.goal.pose.position.y = 15.5;
     }},
  };

  for (const auto& invalid_case : invalid_cases) {
    SCOPED_TRACE(invalid_case.name);
    auto request = makeTestRequest();
    request->k = 1;
    invalid_case.mutate(*request);
    const auto response = callService(executor, client, request);

    ASSERT_NE(response, nullptr);
    EXPECT_FALSE(response->success) << response->message;
    EXPECT_TRUE(response->path_results.empty());
    EXPECT_TRUE(response->debug_nodes.empty());
    EXPECT_EQ(response->result_info.status, PlanningResultInfo::STATUS_INVALID_REQUEST);
    EXPECT_EQ(response->result_info.limits_reached, PlanningResultInfo::LIMIT_NONE);
    EXPECT_EQ(response->result_info.requested_path_count, 1u);
    EXPECT_EQ(response->result_info.found_path_count, 0u);
    EXPECT_EQ(response->result_info.returned_path_count, 0u);
    EXPECT_EQ(response->result_info.expanded_nodes, 0u);
    EXPECT_FALSE(response->result_info.request_satisfied);
    EXPECT_FALSE(response->result_info.search_complete);
    expectStructuredResultInvariants(*response);
    const auto diagnostic = lowercase(response->message);
    EXPECT_NE(diagnostic.find(invalid_case.expected_endpoint), std::string::npos)
      << response->message;
    EXPECT_NE(diagnostic.find(invalid_case.expected_reason), std::string::npos)
      << response->message;
  }

  auto recovery_request = makeTestRequest();
  recovery_request->k = 1;
  const auto recovered = callService(executor, client, recovery_request);
  ASSERT_NE(recovered, nullptr);
  ASSERT_TRUE(recovered->success) << recovered->message;
  ASSERT_EQ(recovered->path_results.size(), 1u);
  EXPECT_EQ(recovered->result_info.status, PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_EQ(recovered->result_info.limits_reached, PlanningResultInfo::LIMIT_NONE);
  EXPECT_EQ(recovered->result_info.requested_path_count, 1u);
  EXPECT_EQ(recovered->result_info.found_path_count, 1u);
  EXPECT_EQ(recovered->result_info.returned_path_count, 1u);
  EXPECT_TRUE(recovered->result_info.request_satisfied);
  EXPECT_TRUE(recovered->result_info.search_complete);
  EXPECT_TRUE(recovered->result_info.output_complete);
  expectStructuredResultInvariants(*recovered);
}

TEST_F(IntegrationTestFixture, RejectedAndNoPathRequestsDoNotPoisonServer) {
  auto node = std::make_shared<rclcpp::Node>("test_invalid_map_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));

  ASSERT_TRUE(client->wait_for_service(2s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto call_service =
    [&](const RaystarService::Request::SharedPtr& request) -> RaystarService::Response::SharedPtr {
    auto future = client->async_send_request(request);
    const auto end_time = std::chrono::steady_clock::now() + 10s;
    while (std::chrono::steady_clock::now() < end_time) {
      if (future.wait_for(100ms) == std::future_status::ready)
        break;
      executor.spin_some(100ms);
    }
    if (future.wait_for(0s) != std::future_status::ready) {
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
  EXPECT_NE(invalid_negative_result->message.find("occupancy value -2"), std::string::npos)
    << invalid_negative_result->message;
  EXPECT_NE(invalid_negative_result->message.find("[-1, 100]"), std::string::npos)
    << invalid_negative_result->message;

  auto invalid_high_request = makeTestRequest();
  invalid_high_request->map.data[0] = static_cast<int8_t>(101);
  auto invalid_high_result = call_service(invalid_high_request);
  ASSERT_NE(invalid_high_result, nullptr);
  EXPECT_FALSE(invalid_high_result->success);
  EXPECT_NE(invalid_high_result->message.find("occupancy value 101"), std::string::npos)
    << invalid_high_result->message;
  EXPECT_NE(invalid_high_result->message.find("[-1, 100]"), std::string::npos)
    << invalid_high_result->message;

  auto blocked_request = makeTestRequest();
  for (unsigned int y = 0; y < blocked_request->map.info.height; ++y)
    blocked_request->map.data[y * blocked_request->map.info.width + 15] = 100;
  auto blocked_result = call_service(blocked_request);
  ASSERT_NE(blocked_result, nullptr);
  EXPECT_FALSE(blocked_result->success);
  EXPECT_NE(blocked_result->message.find("No path exists"), std::string::npos);
  EXPECT_TRUE(blocked_result->debug_nodes.empty());
  EXPECT_EQ(blocked_result->result_info.status, PlanningResultInfo::STATUS_NO_PATH);
  EXPECT_TRUE(blocked_result->result_info.search_complete);
  EXPECT_TRUE(blocked_result->result_info.output_complete);
  EXPECT_FALSE(blocked_result->result_info.request_satisfied);
  EXPECT_EQ(blocked_result->result_info.requested_path_count, 3u);
  EXPECT_EQ(blocked_result->result_info.found_path_count, 0u);
  expectStructuredResultInvariants(*blocked_result);

  auto shared_vertex_request = makeTestRequest();
  std::fill(shared_vertex_request->map.data.begin(), shared_vertex_request->map.data.end(), 0);
  shared_vertex_request->map.data[10 * shared_vertex_request->map.info.width + 10] = 100;
  shared_vertex_request->map.data[11 * shared_vertex_request->map.info.width + 11] = 100;
  auto shared_vertex_result = call_service(shared_vertex_request);
  ASSERT_NE(shared_vertex_result, nullptr);
  EXPECT_FALSE(shared_vertex_result->success);
  EXPECT_NE(shared_vertex_result->message.find("Unsupported obstacle topology"), std::string::npos)
    << shared_vertex_result->message;
  EXPECT_NE(shared_vertex_result->message.find("(11, 11)"), std::string::npos)
    << shared_vertex_result->message;
  EXPECT_TRUE(shared_vertex_result->debug_nodes.empty());
  EXPECT_EQ(shared_vertex_result->result_info.status, PlanningResultInfo::STATUS_FAILED);

  auto valid_result = call_service(makeTestRequest());
  ASSERT_NE(valid_result, nullptr);
  EXPECT_TRUE(valid_result->success);
  EXPECT_GE(valid_result->path_results.size(), 1u);
  expectStructuredResultInvariants(*valid_result);
}

TEST_F(IntegrationTestFixture, RejectsInvalidFramesMetricMetadataAndPosesThenRecovers) {
  auto node = std::make_shared<rclcpp::Node>("test_geometry_contract_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(2s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  using Mutator = std::function<void(RaystarService::Request&)>;
  struct InvalidCase {
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
    {"empty map frame", "frame", [](auto& request) { request.map.header.frame_id.clear(); }},
    {"empty start frame", "frame", [](auto& request) { request.start.header.frame_id.clear(); }},
    {"empty goal frame", "frame", [](auto& request) { request.goal.header.frame_id.clear(); }},
    {"mismatched start frame",
     "frame",
     [](auto& request) { request.start.header.frame_id = "odom"; }},
    {"mismatched goal frame",
     "frame",
     [](auto& request) { request.goal.header.frame_id = "odom"; }},
    {"oversized map frame",
     "256",
     [](auto& request) { request.map.header.frame_id.assign(257, 'm'); }},
    {"zero resolution", "resolution", [](auto& request) { request.map.info.resolution = 0.0f; }},
    {"negative resolution",
     "resolution",
     [](auto& request) { request.map.info.resolution = -1.0f; }},
    {"NaN resolution",
     "resolution",
     [float_nan](auto& request) { request.map.info.resolution = float_nan; }},
    {"infinite resolution",
     "resolution",
     [float_infinity](auto& request) { request.map.info.resolution = float_infinity; }},
    {"non-metric-faithful world transform",
     "metric fidelity",
     [](auto& request) {
       request.map.info.resolution = 10000.0f;
       request.map.info.origin.position.x = 1.0e20;
     }},
    {"tiny resolution collapses next to a large origin",
     "metric fidelity",
     [](auto& request) {
       request.map.info.resolution = 1.0e-12f;
       request.map.info.origin.position.x = 1.0e5;
     }},
    {"origin ULP distorts adjacent cells",
     "metric fidelity",
     [](auto& request) {
       request.map.info.resolution = 0.05f;
       request.map.info.origin.position.x = 1.0e8;
     }},
    {"NaN origin x", "origin", [nan](auto& request) { request.map.info.origin.position.x = nan; }},
    {"infinite origin y",
     "origin",
     [infinity](auto& request) { request.map.info.origin.position.y = infinity; }},
    {"nonzero origin z", "origin", [](auto& request) { request.map.info.origin.position.z = 1.0; }},
    {"NaN origin z", "origin", [nan](auto& request) { request.map.info.origin.position.z = nan; }},
    {"zero origin quaternion",
     "orientation",
     [](auto& request) {
       request.map.info.origin.orientation.x = 0.0;
       request.map.info.origin.orientation.y = 0.0;
       request.map.info.origin.orientation.z = 0.0;
       request.map.info.origin.orientation.w = 0.0;
     }},
    {"non-unit origin quaternion",
     "orientation",
     [](auto& request) { request.map.info.origin.orientation.w = 2.0; }},
    {"nonzero map yaw",
     "orientation",
     [sqrt_half](auto& request) {
       request.map.info.origin.orientation.z = sqrt_half;
       request.map.info.origin.orientation.w = sqrt_half;
     }},
    {"nonzero map roll",
     "orientation",
     [sqrt_half](auto& request) {
       request.map.info.origin.orientation.x = sqrt_half;
       request.map.info.origin.orientation.w = sqrt_half;
     }},
    {"NaN map quaternion",
     "orientation",
     [nan](auto& request) { request.map.info.origin.orientation.z = nan; }},
    {"NaN start x", "start", [nan](auto& request) { request.start.pose.position.x = nan; }},
    {"infinite start y",
     "start",
     [infinity](auto& request) { request.start.pose.position.y = infinity; }},
    {"nonzero start z", "start", [](auto& request) { request.start.pose.position.z = 1.0; }},
    {"infinite goal x",
     "goal",
     [infinity](auto& request) { request.goal.pose.position.x = -infinity; }},
    {"NaN goal y", "goal", [nan](auto& request) { request.goal.pose.position.y = nan; }},
    {"nonzero goal z", "goal", [](auto& request) { request.goal.pose.position.z = -1.0; }},
  };

  for (const auto& invalid_case : invalid_cases) {
    SCOPED_TRACE(invalid_case.name);
    auto request = makeTestRequest();
    invalid_case.mutate(*request);

    auto response = callService(executor, client, request);
    ASSERT_NE(response, nullptr);
    EXPECT_FALSE(response->success) << response->message;
    EXPECT_TRUE(response->path_results.empty());
    EXPECT_TRUE(response->debug_nodes.empty());
    EXPECT_EQ(response->result_info.status, PlanningResultInfo::STATUS_INVALID_REQUEST);
    expectStructuredResultInvariants(*response);
    EXPECT_NE(lowercase(response->message).find(invalid_case.expected_keyword), std::string::npos)
      << response->message;
  }

  auto valid_response = callService(executor, client, makeTestRequest());
  ASSERT_NE(valid_response, nullptr);
  EXPECT_TRUE(valid_response->success) << valid_response->message;
  EXPECT_FALSE(valid_response->path_results.empty());
}

TEST_F(IntegrationTestFixture, RejectsKAboveConfiguredMaximumThenRecovers) {
  auto node = std::make_shared<rclcpp::Node>("test_max_k_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
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
  EXPECT_EQ(excessive_response->result_info.status, PlanningResultInfo::STATUS_INVALID_REQUEST);
  EXPECT_NE(excessive_response->message.find("max_k=100"), std::string::npos)
    << excessive_response->message;

  auto valid_response = callService(executor, client, makeTestRequest());
  ASSERT_NE(valid_response, nullptr);
  EXPECT_TRUE(valid_response->success) << valid_response->message;
  EXPECT_FALSE(valid_response->path_results.empty());
}

TEST_F(IntegrationTestFixture, DynamicParametersAreDescribedAndUpdatedAtomically) {
  const auto spawned = spawnRaystarNode("dynamic_parameters");
  ASSERT_GT(spawned.pid, 0) << "Unable to exec dynamic-parameters raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  ChildProcessGuard parameter_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_dynamic_parameter_client");
  auto service_client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  auto parameter_client =
    std::make_shared<rclcpp::AsyncParametersClient>(node, raystarEndpoint(spawned.node_namespace));
  ASSERT_TRUE(service_client->wait_for_service(5s));
  ASSERT_TRUE(parameter_client->wait_for_service(5s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  struct ExpectedDescriptor {
    const char* name;
    int64_t minimum;
    int64_t maximum;
    bool read_only;
  };
  const int64_t max_size_t_parameter = []() constexpr {
    if constexpr (sizeof(size_t) < sizeof(int64_t))
      return static_cast<int64_t>(std::numeric_limits<size_t>::max());
    return std::numeric_limits<int64_t>::max();
  }
  ();
  const int64_t max_int = static_cast<int64_t>(std::numeric_limits<int>::max());
  const int64_t max_timeout = std::chrono::milliseconds::max().count() - 1;
  const int64_t max_timer_period = std::chrono::duration_cast<std::chrono::milliseconds>(
                                     std::chrono::nanoseconds::max() - std::chrono::milliseconds(1))
                                     .count();
  const std::vector<ExpectedDescriptor> expected_descriptors = {
    {"occupied_threshold", 1, 100, true},
    {"max_k", 1, max_int, false},
    {"max_cost_bounded_paths", 1, max_int, false},
    {"max_multi_goal_count", 1, max_int, false},
    {"max_transition_configurations", 1, max_int, false},
    {"max_transition_pairs", 1, max_int, false},
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
  for (const auto& expected : expected_descriptors) parameter_names.emplace_back(expected.name);

  auto descriptor_future = parameter_client->describe_parameters(parameter_names);
  ASSERT_TRUE(waitForFuture(executor, descriptor_future, 5s));
  const auto descriptors = descriptor_future.get();
  ASSERT_EQ(descriptors.size(), expected_descriptors.size());
  for (size_t i = 0; i < descriptors.size(); ++i) {
    SCOPED_TRACE(expected_descriptors[i].name);
    const auto& descriptor = descriptors[i];
    EXPECT_EQ(descriptor.name, expected_descriptors[i].name);
    EXPECT_EQ(descriptor.type, static_cast<uint8_t>(rclcpp::ParameterType::PARAMETER_INTEGER));
    EXPECT_FALSE(descriptor.description.empty());
    EXPECT_EQ(descriptor.read_only, expected_descriptors[i].read_only);
    ASSERT_EQ(descriptor.integer_range.size(), 1u);
    EXPECT_EQ(descriptor.integer_range.front().from_value, expected_descriptors[i].minimum);
    EXPECT_EQ(descriptor.integer_range.front().to_value, expected_descriptors[i].maximum);
    EXPECT_EQ(descriptor.integer_range.front().step, 1u);
  }

  const std::vector<std::string> initial_parameter_names = {
    "max_k",
    "max_cost_bounded_paths",
    "max_multi_goal_count",
    "max_transition_configurations",
    "max_transition_pairs",
    "occupied_threshold",
    "path_visualization_republish_period_ms",
  };
  auto initial_future = parameter_client->get_parameters(initial_parameter_names);
  ASSERT_TRUE(waitForFuture(executor, initial_future, 5s));
  const auto initial_values = initial_future.get();
  ASSERT_EQ(initial_values.size(), 7u);
  EXPECT_EQ(initial_values[0].as_int(), 100);
  EXPECT_EQ(initial_values[1].as_int(), 1000);
  EXPECT_EQ(initial_values[2].as_int(), 128);
  EXPECT_EQ(initial_values[3].as_int(), 4096);
  EXPECT_EQ(initial_values[4].as_int(), 1000);
  EXPECT_EQ(initial_values[5].as_int(), 99);
  EXPECT_EQ(initial_values[6].as_int(), 2000);

  auto valid_update_future = parameter_client->set_parameters_atomically({
    rclcpp::Parameter("max_k", 1),
  });
  ASSERT_TRUE(waitForFuture(executor, valid_update_future, 5s));
  const auto valid_update = valid_update_future.get();
  ASSERT_TRUE(valid_update.successful) << valid_update.reason;

  auto updated_values_future = parameter_client->get_parameters({"max_k", "occupied_threshold"});
  ASSERT_TRUE(waitForFuture(executor, updated_values_future, 5s));
  const auto updated_values = updated_values_future.get();
  ASSERT_EQ(updated_values.size(), 2u);
  EXPECT_EQ(updated_values[0].as_int(), 1);
  EXPECT_EQ(updated_values[1].as_int(), 99);

  auto excessive_k_request = makeTestRequest();
  excessive_k_request->k = 2;
  auto excessive_k_response = callService(executor, service_client, excessive_k_request);
  ASSERT_NE(excessive_k_response, nullptr);
  EXPECT_FALSE(excessive_k_response->success);
  EXPECT_NE(excessive_k_response->message.find("max_k=1"), std::string::npos)
    << excessive_k_response->message;

  auto below_threshold_response =
    callService(executor, service_client, makeVerticalBarrierRequest(98));
  ASSERT_NE(below_threshold_response, nullptr);
  EXPECT_TRUE(below_threshold_response->success) << below_threshold_response->message;
  auto at_threshold_response =
    callService(executor, service_client, makeVerticalBarrierRequest(99));
  ASSERT_NE(at_threshold_response, nullptr);
  EXPECT_FALSE(at_threshold_response->success);
  EXPECT_EQ(at_threshold_response->result_info.status, PlanningResultInfo::STATUS_NO_PATH);

  auto occupancy_policy_update_future = parameter_client->set_parameters_atomically({
    rclcpp::Parameter("occupied_threshold", 50),
  });
  ASSERT_TRUE(waitForFuture(executor, occupancy_policy_update_future, 5s));
  const auto occupancy_policy_update = occupancy_policy_update_future.get();
  EXPECT_FALSE(occupancy_policy_update.successful);
  EXPECT_FALSE(occupancy_policy_update.reason.empty());
  auto threshold_after_update_future = parameter_client->get_parameters({"occupied_threshold"});
  ASSERT_TRUE(waitForFuture(executor, threshold_after_update_future, 5s));
  const auto threshold_after_update = threshold_after_update_future.get();
  ASSERT_EQ(threshold_after_update.size(), 1u);
  EXPECT_EQ(threshold_after_update.front().as_int(), 99);

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
    rclcpp::Parameter("occupied_threshold", 50),
  });
  ASSERT_TRUE(waitForFuture(executor, mixed_update_future, 5s));
  const auto mixed_update = mixed_update_future.get();
  EXPECT_FALSE(mixed_update.successful);
  EXPECT_FALSE(mixed_update.reason.empty());

  auto unchanged_values_future = parameter_client->get_parameters({"max_k", "occupied_threshold"});
  ASSERT_TRUE(waitForFuture(executor, unchanged_values_future, 5s));
  const auto unchanged_values = unchanged_values_future.get();
  ASSERT_EQ(unchanged_values.size(), 2u);
  EXPECT_EQ(unchanged_values[0].as_int(), 1);
  EXPECT_EQ(unchanged_values[1].as_int(), 99);

  auto visualization_update_future = parameter_client->set_parameters_atomically({
    rclcpp::Parameter("path_visualization_republish_period_ms", 0),
  });
  ASSERT_TRUE(waitForFuture(executor, visualization_update_future, 5s));
  const auto visualization_update = visualization_update_future.get();
  EXPECT_FALSE(visualization_update.successful);
  EXPECT_FALSE(visualization_update.reason.empty());
  auto period_future = parameter_client->get_parameters({"path_visualization_republish_period_ms"});
  ASSERT_TRUE(waitForFuture(executor, period_future, 5s));
  const auto period = period_future.get();
  ASSERT_EQ(period.size(), 1u);
  EXPECT_EQ(period.front().as_int(), 2000);

  auto recovery_update_future = parameter_client->set_parameters_atomically({
    rclcpp::Parameter("max_k", 100),
  });
  ASSERT_TRUE(waitForFuture(executor, recovery_update_future, 5s));
  const auto recovery_update = recovery_update_future.get();
  ASSERT_TRUE(recovery_update.successful) << recovery_update.reason;
  auto recovered_values_future = parameter_client->get_parameters({"max_k", "occupied_threshold"});
  ASSERT_TRUE(waitForFuture(executor, recovered_values_future, 5s));
  const auto recovered_values = recovered_values_future.get();
  ASSERT_EQ(recovered_values.size(), 2u);
  EXPECT_EQ(recovered_values[0].as_int(), 100);
  EXPECT_EQ(recovered_values[1].as_int(), 99);
  auto recovery_response = callService(executor, service_client, makeTestRequest());
  ASSERT_NE(recovery_response, nullptr);
  EXPECT_TRUE(recovery_response->success) << recovery_response->message;
  expectStructuredResultInvariants(*recovery_response);
}

TEST_F(IntegrationTestFixture, ConfiguredOccupancyThresholdIsApplied) {
  const auto spawned = spawnRaystarNode("occupancy_threshold_50", {"occupied_threshold:=50"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec occupancy-threshold raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  ChildProcessGuard threshold_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_configured_occupancy_threshold_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  auto parameter_client =
    std::make_shared<rclcpp::AsyncParametersClient>(node, raystarEndpoint(spawned.node_namespace));
  ASSERT_TRUE(client->wait_for_service(5s));
  ASSERT_TRUE(parameter_client->wait_for_service(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto configured_threshold_future = parameter_client->get_parameters({"occupied_threshold"});
  ASSERT_TRUE(waitForFuture(executor, configured_threshold_future, 5s));
  const auto configured_threshold = configured_threshold_future.get();
  ASSERT_EQ(configured_threshold.size(), 1u);
  EXPECT_EQ(configured_threshold.front().as_int(), 50);

  auto below_threshold = callService(executor, client, makeVerticalBarrierRequest(49));
  ASSERT_NE(below_threshold, nullptr);
  EXPECT_TRUE(below_threshold->success) << below_threshold->message;

  auto at_threshold = callService(executor, client, makeVerticalBarrierRequest(50));
  ASSERT_NE(at_threshold, nullptr);
  EXPECT_FALSE(at_threshold->success);
  EXPECT_NE(at_threshold->message.find("No path exists"), std::string::npos)
    << at_threshold->message;
}

TEST_F(IntegrationTestFixture, RejectsInvalidOccupancyThresholdAtStartup) {
  struct InvalidThresholdCase {
    const char* namespace_remap;
    const char* parameter_override;
  };
  const std::vector<InvalidThresholdCase> invalid_cases = {
    {"occupancy_threshold_zero", "occupied_threshold:=0"},
    {"occupancy_threshold_101", "occupied_threshold:=101"},
  };

  for (const auto& invalid_case : invalid_cases) {
    SCOPED_TRACE(invalid_case.parameter_override);
    const auto spawned =
      spawnRaystarNode(invalid_case.namespace_remap, {invalid_case.parameter_override});
    ASSERT_GT(spawned.pid, 0) << "Unable to exec invalid-parameter raystar_node: "
                              << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                       : "fork/pipe failure");
    ChildProcessGuard invalid_process(spawned.pid);
    int status = 0;
    ASSERT_TRUE(waitForChildExit(spawned.pid, status, 5s))
      << "raystar_node accepted " << invalid_case.parameter_override;
    invalid_process.disarm();
    ASSERT_TRUE(WIFEXITED(status));
    EXPECT_EQ(WEXITSTATUS(status), 1);
  }
}

TEST_F(IntegrationTestFixture, MaxNodesStopsExpansionAndLimitedServerRecovers) {
  const auto spawned = spawnRaystarNode(
    "resource_limited", {"max_nodes:=1", "planning_timeout_ms:=60000", "max_debug_nodes:=1"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec max-nodes raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
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
  EXPECT_EQ(limited_response->result_info.status, PlanningResultInfo::STATUS_PARTIAL_SEARCH);
  EXPECT_NE(limited_response->result_info.limits_reached & PlanningResultInfo::LIMIT_MAX_NODES, 0u);
  EXPECT_FALSE(limited_response->result_info.search_complete);
  EXPECT_TRUE(limited_response->result_info.output_complete);
  EXPECT_FALSE(limited_response->result_info.request_satisfied);
  EXPECT_EQ(limited_response->result_info.requested_path_count, 3u);
  EXPECT_EQ(limited_response->result_info.found_path_count, 0u);
  EXPECT_EQ(limited_response->result_info.returned_path_count, 0u);
  expectStructuredResultInvariants(*limited_response);
  EXPECT_NE(limited_response->message.find("max_nodes=1"), std::string::npos)
    << limited_response->message;

  auto bounded_limited_request = makeTestRequest();
  bounded_limited_request->search_mode = RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH;
  bounded_limited_request->k = 0;
  bounded_limited_request->max_path_length = 1000.0;
  const auto bounded_limited_response = callService(executor, client, bounded_limited_request);
  ASSERT_NE(bounded_limited_response, nullptr);
  EXPECT_FALSE(bounded_limited_response->success);
  EXPECT_EQ(bounded_limited_response->result_info.status,
            PlanningResultInfo::STATUS_PARTIAL_SEARCH);
  EXPECT_NE(
    bounded_limited_response->result_info.limits_reached & PlanningResultInfo::LIMIT_MAX_NODES, 0u);
  EXPECT_FALSE(bounded_limited_response->result_info.search_complete);
  EXPECT_FALSE(bounded_limited_response->result_info.cost_bound_exhausted);
  EXPECT_TRUE(bounded_limited_response->result_info.output_complete);
  EXPECT_FALSE(bounded_limited_response->result_info.request_satisfied);
  expectStructuredResultInvariants(*bounded_limited_response);

  auto root_solution_request = makeTestRequest();
  std::fill(root_solution_request->map.data.begin(), root_solution_request->map.data.end(), 0);
  root_solution_request->k = 1;
  root_solution_request->include_debug = true;
  auto recovered_response = callService(executor, client, root_solution_request);
  ASSERT_NE(recovered_response, nullptr);
  EXPECT_TRUE(recovered_response->success) << recovered_response->message;
  EXPECT_EQ(recovered_response->path_results.size(), 1u);
  EXPECT_EQ(recovered_response->debug_nodes.size(), 1u);
  EXPECT_EQ(recovered_response->result_info.status, PlanningResultInfo::STATUS_COMPLETE);
  EXPECT_TRUE(recovered_response->result_info.request_satisfied);
  expectStructuredResultInvariants(*recovered_response);
}

TEST_F(IntegrationTestFixture, MaxCostBoundedPathsReturnsUncertifiedPartialSearch) {
  const auto spawned = spawnRaystarNode("bounded_path_limited", {"max_cost_bounded_paths:=1"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec max-paths raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  ChildProcessGuard limited_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_max_cost_bounded_paths_client");
  auto client = node->create_client<RaystarService>(
    raystarEndpoint(spawned.node_namespace, "get_raystar_paths"));
  ASSERT_TRUE(client->wait_for_service(5s));
  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);

  auto request = makeTestRequest();
  request->search_mode = RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH;
  request->k = 0;
  request->max_path_length = 100.0;
  const auto response = callService(executor, client, request);

  ASSERT_NE(response, nullptr);
  ASSERT_TRUE(response->success) << response->message;
  ASSERT_EQ(response->path_results.size(), 1u);
  EXPECT_EQ(response->result_info.status, PlanningResultInfo::STATUS_PARTIAL_SEARCH);
  EXPECT_NE(response->result_info.limits_reached & PlanningResultInfo::LIMIT_MAX_PATHS, 0u);
  EXPECT_FALSE(response->result_info.cost_bound_exhausted);
  EXPECT_FALSE(response->result_info.search_complete);
  EXPECT_TRUE(response->result_info.output_complete);
  EXPECT_FALSE(response->result_info.request_satisfied);
  EXPECT_EQ(response->result_info.requested_path_count, 0u);
  EXPECT_EQ(response->result_info.found_path_count, 1u);
  EXPECT_EQ(response->result_info.returned_path_count, 1u);
  EXPECT_LE(response->path_results.front().cost, request->max_path_length);
  expectStructuredResultInvariants(*response);
  EXPECT_NE(response->message.find("max_cost_bounded_paths=1"), std::string::npos)
    << response->message;
}

TEST_F(IntegrationTestFixture, PlanningTimeoutHasStructuredLimitStatus) {
  const auto spawned = spawnRaystarNode("planning_timeout", {"planning_timeout_ms:=1"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec planning-timeout raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
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
  EXPECT_EQ(response->result_info.status, PlanningResultInfo::STATUS_PARTIAL_SEARCH);
  EXPECT_NE(response->result_info.limits_reached & PlanningResultInfo::LIMIT_TIMEOUT, 0u);
  EXPECT_FALSE(response->result_info.search_complete);
  EXPECT_FALSE(response->result_info.request_satisfied);
  EXPECT_TRUE(response->path_results.empty());
  expectStructuredResultInvariants(*response);
}

TEST_F(IntegrationTestFixture, MapCellAdmissionHappensBeforePlannerAndRecovers) {
  const auto spawned = spawnRaystarNode("map_cell_limited", {"max_map_cells:=899"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec map-cell-limit raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
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
  EXPECT_NE(rejected->message.find("max_map_cells=899"), std::string::npos) << rejected->message;

  auto smaller = makeTestRequest();
  smaller->map.info.width = 29;
  smaller->map.data.assign(29 * 30, 0);
  for (unsigned int y = 10; y < 20; ++y)
    for (unsigned int x = 10; x < 20; ++x) smaller->map.data[y * 29 + x] = 100;
  auto recovered = callService(executor, client, smaller);
  ASSERT_NE(recovered, nullptr);
  EXPECT_TRUE(recovered->success) << recovered->message;
  EXPECT_FALSE(recovered->path_results.empty());
}

TEST_F(IntegrationTestFixture, MapByteAdmissionUsesWorkingSetEstimate) {
  const auto spawned =
    spawnRaystarNode("map_byte_limited", {"max_map_cells:=1000", "max_map_bytes:=28000"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec map-byte-limit raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
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
  EXPECT_NE(rejected->message.find("max_map_bytes=28000"), std::string::npos) << rejected->message;

  auto smaller = makeTestRequest();  // 29 * 30 * 32 = 27840 bytes
  smaller->map.info.width = 29;
  smaller->map.data.assign(29 * 30, 0);
  for (unsigned int y = 10; y < 20; ++y)
    for (unsigned int x = 10; x < 20; ++x) smaller->map.data[y * 29 + x] = 100;
  auto recovered = callService(executor, client, smaller);
  ASSERT_NE(recovered, nullptr);
  EXPECT_TRUE(recovered->success) << recovered->message;
}

TEST_F(IntegrationTestFixture, PathInterpolationLimitRejectsOversizedPath) {
  const auto spawned = spawnRaystarNode("path_points_limited", {"max_path_points:=10"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec path-point-limit raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
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
  EXPECT_EQ(rejected->result_info.status, PlanningResultInfo::STATUS_PARTIAL_OUTPUT);
  EXPECT_NE(rejected->result_info.limits_reached & PlanningResultInfo::LIMIT_MAX_PATH_POINTS, 0u);
  EXPECT_TRUE(rejected->result_info.search_complete);
  EXPECT_FALSE(rejected->result_info.output_complete);
  EXPECT_EQ(rejected->result_info.requested_path_count, 1u);
  EXPECT_EQ(rejected->result_info.found_path_count, 1u);
  EXPECT_EQ(rejected->result_info.returned_path_count, 0u);
  expectStructuredResultInvariants(*rejected);
  EXPECT_NE(rejected->message.find("max_path_points=10"), std::string::npos) << rejected->message;

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

TEST_F(IntegrationTestFixture, PathPointLimitIsCumulativeAcrossResponse) {
  const auto spawned = spawnRaystarNode("total_path_points_limited", {"max_path_points:=30"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec cumulative path-point-limit raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  ChildProcessGuard limited_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_total_path_points_limit_client");
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
  EXPECT_EQ(response->result_info.status, PlanningResultInfo::STATUS_PARTIAL_OUTPUT);
  EXPECT_NE(response->result_info.limits_reached & PlanningResultInfo::LIMIT_MAX_PATH_POINTS, 0u);
  EXPECT_TRUE(response->result_info.search_complete);
  EXPECT_FALSE(response->result_info.output_complete);
  EXPECT_FALSE(response->result_info.request_satisfied);
  EXPECT_EQ(response->result_info.requested_path_count, 2u);
  EXPECT_EQ(response->result_info.found_path_count, 2u);
  EXPECT_EQ(response->result_info.returned_path_count, 1u);
  expectStructuredResultInvariants(*response);
  EXPECT_NE(response->message.find("per-response max_path_points=30"), std::string::npos)
    << response->message;

  auto bounded_request = makeTestRequest();
  bounded_request->search_mode = RaystarService::Request::SEARCH_MODE_ALL_WITHIN_LENGTH;
  bounded_request->k = 0;
  bounded_request->max_path_length = 100.0;
  const auto bounded_response = callService(executor, client, bounded_request);
  ASSERT_NE(bounded_response, nullptr);
  ASSERT_TRUE(bounded_response->success) << bounded_response->message;
  ASSERT_EQ(bounded_response->path_results.size(), 1u) << bounded_response->message;
  EXPECT_EQ(bounded_response->result_info.status, PlanningResultInfo::STATUS_PARTIAL_OUTPUT);
  EXPECT_TRUE(bounded_response->result_info.cost_bound_exhausted);
  EXPECT_TRUE(bounded_response->result_info.search_complete);
  EXPECT_FALSE(bounded_response->result_info.output_complete);
  EXPECT_FALSE(bounded_response->result_info.request_satisfied);
  EXPECT_EQ(bounded_response->result_info.found_path_count, 2u);
  EXPECT_EQ(bounded_response->result_info.returned_path_count, 1u);
  expectStructuredResultInvariants(*bounded_response);
}

TEST_F(IntegrationTestFixture, DebugNodeLimitBoundsResponseArrays) {
  const auto spawned = spawnRaystarNode("debug_nodes_limited", {"max_debug_nodes:=0"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec debug-node-limit raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
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
  EXPECT_NE(response->message.find("Debug output limited"), std::string::npos) << response->message;
}

TEST_F(IntegrationTestFixture, ResponseByteLimitBoundsPathsAndRecovers) {
  const auto spawned = spawnRaystarNode("response_bytes_limited", {"max_response_bytes:=1024"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec response-byte-limit raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
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
  EXPECT_EQ(rejected->result_info.status, PlanningResultInfo::STATUS_PARTIAL_OUTPUT);
  EXPECT_NE(rejected->result_info.limits_reached & PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES,
            0u);
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

TEST_F(IntegrationTestFixture, ReservedNoDeadlineTimeoutFailsAtStartup) {
  const auto spawned =
    spawnRaystarNode("invalid_timeout", {"planning_timeout_ms:=9223372036854775807"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec invalid-timeout raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  ChildProcessGuard invalid_timeout_process(spawned.pid);
  int status = 0;
  ASSERT_TRUE(waitForChildExit(spawned.pid, status, 5s));
  invalid_timeout_process.disarm();
  ASSERT_TRUE(WIFEXITED(status));
  EXPECT_EQ(WEXITSTATUS(status), 1);
}

TEST_F(IntegrationTestFixture, MarkerSnapshotsReplaceClearAndRemainDurable) {
  const auto spawned = spawnRaystarNode("marker_lifecycle",
                                        {"planning_timeout_ms:=60000",
                                         "path_visualization_republish_period_ms:=200",
                                         "max_debug_nodes:=100"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec marker-lifecycle raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
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
  for (size_t path_index = 0; path_index < first_response->path_results.size(); ++path_index) {
    SCOPED_TRACE("returned path " + std::to_string(path_index));
    expectIndependentCollisionFreePath(first_request->map,
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
    raystarEndpoint(spawned.node_namespace, "non_homotopic_paths"),
    marker_qos,
    [&](MarkerArray::ConstSharedPtr message) { path_messages.push_back(std::move(message)); }));
  subscriptions.push_back(node->create_subscription<MarkerArray>(
    raystarEndpoint(spawned.node_namespace, "poly_obstacles"),
    marker_qos,
    [&](MarkerArray::ConstSharedPtr message) { obstacle_messages.push_back(std::move(message)); }));
  subscriptions.push_back(node->create_subscription<MarkerArray>(
    raystarEndpoint(spawned.node_namespace, "cdt"),
    marker_qos,
    [&](MarkerArray::ConstSharedPtr message) { cdt_messages.push_back(std::move(message)); }));
  subscriptions.push_back(node->create_subscription<MarkerArray>(
    raystarEndpoint(spawned.node_namespace, "debug_tree"),
    marker_qos,
    [&](MarkerArray::ConstSharedPtr message) { tree_messages.push_back(std::move(message)); }));

  const auto marker_stamp_ns = [](const MarkerArray& array) {
    const auto& stamp = array.markers.front().header.stamp;
    return static_cast<int64_t>(stamp.sec) * 1000000000LL + static_cast<int64_t>(stamp.nanosec);
  };
  const auto is_complete_replacement = [](const MarkerArray& array) {
    if (array.markers.size() < 2 ||
        array.markers.front().action != visualization_msgs::msg::Marker::DELETEALL) {
      return false;
    }
    return std::all_of(array.markers.begin() + 1, array.markers.end(), [](const auto& marker) {
      return marker.action == visualization_msgs::msg::Marker::ADD;
    });
  };
  const auto is_replacement_snapshot =
    [&](const MarkerArray& array, int64_t min_stamp_ns, int64_t max_stamp_ns) {
      if (!is_complete_replacement(array))
        return false;
      const int64_t stamp_ns = marker_stamp_ns(array);
      if (stamp_ns < min_stamp_ns || stamp_ns > max_stamp_ns)
        return false;
      return true;
    };
  const auto is_clear_snapshot =
    [&](const MarkerArray& array, int64_t min_stamp_ns, int64_t max_stamp_ns) {
      return array.markers.size() == 1 &&
             array.markers.front().action == visualization_msgs::msg::Marker::DELETEALL &&
             marker_stamp_ns(array) >= min_stamp_ns && marker_stamp_ns(array) <= max_stamp_ns;
    };
  const auto find_message = [](const MarkerHistory& history,
                               const auto& predicate) -> MarkerArray::ConstSharedPtr {
    const auto found = std::find_if(
      history.begin(), history.end(), [&](const auto& message) { return predicate(*message); });
    return found == history.end() ? nullptr : *found;
  };
  const auto spin_until = [&](const auto& predicate, std::chrono::milliseconds timeout) {
    const auto end_time = std::chrono::steady_clock::now() + timeout;
    while (std::chrono::steady_clock::now() < end_time) {
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
  ASSERT_TRUE(spin_until(
    [&]() {
      retained_paths = find_message(path_messages, [&](const auto& array) {
        return is_complete_replacement(array) &&
               array.markers.size() == first_response->path_results.size() + 1;
      });
      retained_obstacles = find_message(obstacle_messages, [&](const auto& array) {
        return is_replacement_snapshot(array, first_request_start_ns, subscription_start_ns);
      });
      retained_cdt = find_message(cdt_messages, [&](const auto& array) {
        return is_replacement_snapshot(array, first_request_start_ns, subscription_start_ns);
      });
      retained_tree = find_message(tree_messages, [&](const auto& array) {
        return is_replacement_snapshot(array, first_request_start_ns, subscription_start_ns);
      });
      return retained_paths && retained_obstacles && retained_cdt && retained_tree;
    },
    3s))
    << "Late subscribers did not receive retained marker snapshots";

  for (const auto* snapshot :
       {retained_paths.get(), retained_obstacles.get(), retained_cdt.get(), retained_tree.get()}) {
    ASSERT_NE(snapshot, nullptr);
    for (const auto& marker : snapshot->markers) EXPECT_EQ(marker.header.frame_id, "map");
  }

  // RViz deletes a disabled path namespace locally and needs a fresh message
  // before it can display that namespace again.  Only the already-built path
  // snapshot should be repeated; the much larger geometry/debug snapshots
  // must remain one-shot transient-local publications.
  path_messages.clear();
  obstacle_messages.clear();
  cdt_messages.clear();
  tree_messages.clear();
  ASSERT_TRUE(spin_until([&]() { return path_messages.size() >= 2; }, 1500ms))
    << "Cached path visualization was not periodically republished";
  EXPECT_TRUE(obstacle_messages.empty());
  EXPECT_TRUE(cdt_messages.empty());
  EXPECT_TRUE(tree_messages.empty());
  for (const auto& repeated_paths : path_messages) {
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
  ASSERT_TRUE(spin_until(
    [&]() {
      fewer_paths = find_message(path_messages, [&](const auto& array) {
        return is_complete_replacement(array) && array.markers.size() == 2;
      });
      return fewer_paths != nullptr;
    },
    2s));
  ASSERT_EQ(fewer_paths->markers.size(), 2u);
  EXPECT_EQ(fewer_paths->markers[1].ns, "path_1");
  EXPECT_EQ(fewer_paths->markers[1].id, 1);

  // The periodic cache must be replaced transactionally with the new K=1
  // result rather than continuing to replay the previous K=3 snapshot.
  path_messages.clear();
  ASSERT_TRUE(spin_until([&]() { return !path_messages.empty(); }, 1s));
  for (const auto& repeated_paths : path_messages) {
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
             return is_clear_snapshot(array, invalid_request_start_ns, invalid_response_end_ns);
           }) != nullptr;
  };
  ASSERT_TRUE(spin_until(
    [&]() {
      return has_clear_for_request(path_messages) && has_clear_for_request(obstacle_messages) &&
             has_clear_for_request(cdt_messages) && has_clear_for_request(tree_messages);
    },
    2s));

  // Observe by receive order, not by Marker header stamp: cached snapshots
  // deliberately keep their original stamp, so a stale replay after this
  // DELETEALL would otherwise evade a timestamp-based assertion.
  path_messages.clear();
  obstacle_messages.clear();
  cdt_messages.clear();
  tree_messages.clear();
  const auto no_resurrection_end = std::chrono::steady_clock::now() + 800ms;
  while (std::chrono::steady_clock::now() < no_resurrection_end) {
    executor.spin_some(50ms);
    std::this_thread::sleep_for(10ms);
  }
  const auto contains_add_after_failure = [](const MarkerHistory& history) {
    return std::any_of(history.begin(), history.end(), [](const auto& message) {
      return std::any_of(message->markers.begin(), message->markers.end(), [](const auto& marker) {
        return marker.action == visualization_msgs::msg::Marker::ADD;
      });
    });
  };
  EXPECT_FALSE(contains_add_after_failure(path_messages));
  EXPECT_FALSE(contains_add_after_failure(obstacle_messages));
  EXPECT_FALSE(contains_add_after_failure(cdt_messages));
  EXPECT_FALSE(contains_add_after_failure(tree_messages));
}

TEST_F(IntegrationTestFixture, ZeroPathRepublishPeriodKeepsOnlyDurableSnapshot) {
  const auto spawned =
    spawnRaystarNode("marker_refresh_disabled", {"path_visualization_republish_period_ms:=0"});
  ASSERT_GT(spawned.pid, 0) << "Unable to exec marker-refresh raystar_node: "
                            << (spawned.launch_errno ? std::strerror(spawned.launch_errno)
                                                     : "fork/pipe failure");
  ChildProcessGuard marker_process(spawned.pid);

  auto node = std::make_shared<rclcpp::Node>("test_marker_refresh_disabled_client");
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
    [&](MarkerArray::ConstSharedPtr message) { path_messages.push_back(std::move(message)); });
  (void)subscription;

  const auto retained_deadline = std::chrono::steady_clock::now() + 3s;
  while (path_messages.empty() && std::chrono::steady_clock::now() < retained_deadline) {
    executor.spin_some(50ms);
    std::this_thread::sleep_for(10ms);
  }
  ASSERT_FALSE(path_messages.empty()) << "Late subscriber did not receive a retained path snapshot";
  ASSERT_GT(path_messages.back()->markers.size(), 1u);
  EXPECT_EQ(path_messages.back()->markers.front().action,
            visualization_msgs::msg::Marker::DELETEALL);

  // Allow any discovery-time retained delivery already in flight to drain,
  // then require the observed count to stay fixed for longer than the normal
  // two-second default period.  This is robust to an RMW delivering the same
  // retained sample more than once during endpoint matching.
  const auto drain_deadline = std::chrono::steady_clock::now() + 200ms;
  while (std::chrono::steady_clock::now() < drain_deadline) {
    executor.spin_some(50ms);
    std::this_thread::sleep_for(10ms);
  }
  const size_t retained_count = path_messages.size();
  const auto no_refresh_deadline = std::chrono::steady_clock::now() + 2300ms;
  while (std::chrono::steady_clock::now() < no_refresh_deadline) {
    executor.spin_some(50ms);
    std::this_thread::sleep_for(10ms);
  }
  EXPECT_EQ(path_messages.size(), retained_count)
    << "A zero republish period must disable periodic path traffic";
}

TEST_F(IntegrationTestFixture, PathRepublishPeriodDefaultsToTwoSeconds) {
  auto node = std::make_shared<rclcpp::Node>("test_path_refresh_default_parameter_client");
  auto client =
    std::make_shared<rclcpp::AsyncParametersClient>(node, raystarEndpoint(main_namespace_));
  ASSERT_TRUE(client->wait_for_service(5s));

  rclcpp::executors::SingleThreadedExecutor executor;
  executor.add_node(node);
  auto future = client->get_parameters({"path_visualization_republish_period_ms"});
  ASSERT_TRUE(waitForFuture(executor, future, 5s));
  const auto parameters = future.get();
  ASSERT_EQ(parameters.size(), 1u);
  EXPECT_EQ(parameters.front().get_type(), rclcpp::ParameterType::PARAMETER_INTEGER);
  EXPECT_EQ(parameters.front().as_int(), 2000);
}

TEST_F(IntegrationTestFixture, CustomMapFramePropagatesToEveryReturnedPathPose) {
  auto node = std::make_shared<rclcpp::Node>("test_custom_frame_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
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
  for (const auto& path_result : response->path_results) {
    const auto& path = path_result.path;
    EXPECT_EQ(path.header.frame_id, custom_frame);
    ASSERT_FALSE(path.poses.empty());
    for (const auto& pose : path.poses) {
      EXPECT_EQ(pose.header.frame_id, custom_frame);
      EXPECT_TRUE(std::isfinite(pose.pose.position.x));
      EXPECT_TRUE(std::isfinite(pose.pose.position.y));
    }
  }
}

TEST_F(IntegrationTestFixture, NegativeIdentityMapQuaternionIsAccepted) {
  auto node = std::make_shared<rclcpp::Node>("test_equivalent_identity_client");
  auto client =
    node->create_client<RaystarService>(raystarEndpoint(main_namespace_, "get_raystar_paths"));
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

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
