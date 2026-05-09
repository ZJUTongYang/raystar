#include <gtest/gtest.h>
#include <rclcpp/rclcpp.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <raystar_interfaces/srv/get_raystar_paths.hpp>

#include <sys/wait.h>
#include <unistd.h>

#include <chrono>
#include <thread>

using namespace std::chrono_literals;

static nav_msgs::msg::OccupancyGrid makeTestGrid()
{
  nav_msgs::msg::OccupancyGrid grid;
  grid.info.width = 30;
  grid.info.height = 30;
  grid.info.resolution = 1.0f;
  grid.info.origin.position.x = 0.0;
  grid.info.origin.position.y = 0.0;
  grid.data.resize(900, 0);

  for (unsigned int y = 10; y < 20; ++y)
    for (unsigned int x = 10; x < 20; ++x)
      grid.data[y * 30 + x] = 100;

  return grid;
}

class IntegrationTestFixture : public ::testing::Test
{
protected:
  static void SetUpTestSuite()
  {
    rclcpp::init(0, nullptr);

    pid_ = fork();
    if (pid_ == 0) {
      execl(RAYSTAR_NODE_PATH, "raystar_node", nullptr);
      _exit(1);
    }

    auto node = std::make_shared<rclcpp::Node>("test_wait_node");
    auto client = node->create_client<raystar_interfaces::srv::GetRaystarPaths>(
      "/raystar/get_raystar_paths");

    bool found = false;
    for (int i = 0; i < 50; ++i) {
      if (client->wait_for_service(100ms)) { found = true; break; }
    }
    ASSERT_TRUE(found) << "Service /raystar/get_raystar_paths not available after 5s";
  }

  static void TearDownTestSuite()
  {
    if (pid_ > 0) {
      kill(pid_, SIGKILL);
      int status;
      waitpid(pid_, &status, 0);
    }
    auto ctx = rclcpp::contexts::get_global_default_context();
    if (ctx && ctx->is_valid()) {
      ctx->shutdown("test teardown");
    }
  }

  static pid_t pid_;
};

pid_t IntegrationTestFixture::pid_ = 0;

TEST_F(IntegrationTestFixture, ServiceCallReturnsPaths)
{
  auto node = std::make_shared<rclcpp::Node>("test_integration_client");
  auto client = node->create_client<raystar_interfaces::srv::GetRaystarPaths>(
    "/raystar/get_raystar_paths");

  ASSERT_TRUE(client->wait_for_service(2s));

  auto request = std::make_shared<raystar_interfaces::srv::GetRaystarPaths::Request>();
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
  EXPECT_GE(result->paths.size(), 1u);
  EXPECT_EQ(result->paths.size(), result->path_costs.size());
}

int main(int argc, char** argv)
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
