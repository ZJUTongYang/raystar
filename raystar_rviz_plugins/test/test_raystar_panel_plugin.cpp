#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QCoreApplication>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QTableWidget>
#include <QString>
#include <gtest/gtest.h>
#include <geometry_msgs/msg/point_stamped.hpp>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <raystar_interfaces/action/plan_raystar_goal_set.hpp>
#include <raystar_interfaces/action/plan_raystar_paths.hpp>
#include <raystar_interfaces/environment_identity.hpp>
#include <raystar_interfaces/map_identity.hpp>
#include <raystar_interfaces/msg/goal_path_result.hpp>
#include <raystar_interfaces/msg/multi_goal_planning_result_info.hpp>
#include <raystar_interfaces/msg/path_result.hpp>
#include <raystar_interfaces/msg/planning_result_info.hpp>
#include <rclcpp/rclcpp.hpp>
#include <rclcpp_action/rclcpp_action.hpp>
#include <rviz_common/display_context.hpp>
#include <pluginlib/class_loader.hpp>
#include <rviz_common/config.hpp>
#include <rviz_common/panel.hpp>
#include <rviz_common/ros_integration/ros_node_abstraction_iface.hpp>

#include "raystar_rviz_plugins/raystar_panel.h"

namespace {

using PlanningAction = raystar_interfaces::action::PlanRaystarPaths;
using GoalHandle = rclcpp_action::ServerGoalHandle<PlanningAction>;
using GoalSetAction = raystar_interfaces::action::PlanRaystarGoalSet;
using GoalSetGoalHandle = rclcpp_action::ServerGoalHandle<GoalSetAction>;

class TestRosNodeAbstraction final : public rviz_common::ros_integration::RosNodeAbstractionIface {
public:
  explicit TestRosNodeAbstraction(rclcpp::Node::SharedPtr node) : node_(std::move(node)) {}

  std::string get_node_name() const override {
    return node_ ? node_->get_name() : std::string();
  }

  std::map<std::string, std::vector<std::string>> get_topic_names_and_types() const override {
    return node_ ? node_->get_topic_names_and_types()
                 : std::map<std::string, std::vector<std::string>>{};
  }

  rclcpp::Node::SharedPtr get_raw_node() override {
    return node_;
  }

private:
  rclcpp::Node::SharedPtr node_;
};

// RaystarPanel only needs the ROS-node abstraction from DisplayContext.  The
// remaining methods are deliberately inert, which keeps this test independent
// of Ogre/rendering initialization while still exercising Panel::initialize().
class TestDisplayContext final : public rviz_common::DisplayContext {
public:
  explicit TestDisplayContext(std::shared_ptr<TestRosNodeAbstraction> abstraction)
    : abstraction_(std::move(abstraction)) {}

  Ogre::SceneManager* getSceneManager() const override {
    return nullptr;
  }
  rviz_common::WindowManagerInterface* getWindowManager() const override {
    return nullptr;
  }
  std::shared_ptr<rviz_common::interaction::SelectionManagerIface> getSelectionManager()
    const override {
    return nullptr;
  }
  std::shared_ptr<rviz_common::interaction::HandlerManagerIface> getHandlerManager()
    const override {
    return nullptr;
  }
  std::shared_ptr<rviz_common::interaction::ViewPickerIface> getViewPicker() const override {
    return nullptr;
  }
  rviz_common::FrameManagerIface* getFrameManager() const override {
    return nullptr;
  }
  QString getFixedFrame() const override {
    return QStringLiteral("map");
  }
  uint64_t getFrameCount() const override {
    return 0;
  }
  rviz_common::DisplayFactory* getDisplayFactory() const override {
    return nullptr;
  }
  rviz_common::ros_integration::RosNodeAbstractionIface::WeakPtr getRosNodeAbstraction()
    const override {
    return abstraction_;
  }
  void handleChar(QKeyEvent*, rviz_common::RenderPanel*) override {}
  void handleMouseEvent(const rviz_common::ViewportMouseEvent&) override {}
  rviz_common::ToolManager* getToolManager() const override {
    return nullptr;
  }
  rviz_common::ViewManager* getViewManager() const override {
    return nullptr;
  }
  rviz_common::transformation::TransformationManager* getTransformationManager() override {
    return nullptr;
  }
  rviz_common::DisplayGroup* getRootDisplayGroup() const override {
    return nullptr;
  }
  uint32_t getDefaultVisibilityBit() const override {
    return 0;
  }
  rviz_common::BitAllocator* visibilityBits() override {
    return nullptr;
  }
  void setStatus(const QString&) override {}
  QString getHelpPath() const override {
    return {};
  }
  std::shared_ptr<rclcpp::Clock> getClock() override {
    return std::make_shared<rclcpp::Clock>(RCL_STEADY_TIME);
  }
  void lockRender() override {}
  void unlockRender() override {}
  void queueRender() override {}

private:
  std::shared_ptr<TestRosNodeAbstraction> abstraction_;
};

enum class ServerMode { success, mismatched_map, delayed, inner_bound_mismatch, partial_output };

class TestActionServer {
public:
  TestActionServer(std::string action_name,
                   ServerMode mode,
                   std::chrono::milliseconds execution_delay = std::chrono::milliseconds::zero())
    : node_(rclcpp::Node::make_shared(uniqueName("server")))
    , client_node_(rclcpp::Node::make_shared(uniqueName("client")))
    , action_name_(std::move(action_name))
    , mode_(mode)
    , execution_delay_(execution_delay) {
    server_ = rclcpp_action::create_server<PlanningAction>(
      node_,
      action_name_,
      [this](const rclcpp_action::GoalUUID&, std::shared_ptr<const PlanningAction::Goal> goal) {
        {
          std::lock_guard<std::mutex> lock(goal_mutex_);
          last_goal_ = *goal;
        }
        goal_count_.fetch_add(1);
        return rclcpp_action::GoalResponse::ACCEPT_AND_EXECUTE;
      },
      [this](const std::shared_ptr<GoalHandle>&) {
        cancel_count_.fetch_add(1);
        return rclcpp_action::CancelResponse::ACCEPT;
      },
      [this](const std::shared_ptr<GoalHandle>& goal_handle) {
        std::lock_guard<std::mutex> lock(worker_mutex_);
        if (stopping_) {
          return;
        }
        workers_.emplace_back([this, goal_handle]() noexcept {
          try {
            const auto deadline = std::chrono::steady_clock::now() + execution_delay_;
            while (std::chrono::steady_clock::now() < deadline) {
              if (goal_handle->is_canceling()) {
                finishCanceled(goal_handle);
                return;
              }
              std::unique_lock<std::mutex> lock(worker_mutex_);
              if (worker_cv_.wait_for(
                    lock, std::chrono::milliseconds(5), [this]() { return stopping_; })) {
                return;
              }
            }
            {
              std::lock_guard<std::mutex> lock(worker_mutex_);
              if (stopping_) {
                return;
              }
            }
            if (goal_handle->is_canceling()) {
              finishCanceled(goal_handle);
              return;
            }
            auto result = std::make_shared<PlanningAction::Result>();
            result->success = mode_ == ServerMode::success;
            result->result_info.status =
              raystar_interfaces::msg::PlanningResultInfo::STATUS_COMPLETE;
            result->result_info.request_satisfied = true;
            result->result_info.search_complete = true;
            result->result_info.output_complete = true;
            result->result_info.requested_path_count = 1;
            result->result_info.found_path_count = 1;
            result->result_info.returned_path_count = 1;
            result->result_info.map_id = goal_handle->get_goal()->map_id;
            if (mode_ == ServerMode::mismatched_map) {
              result->result_info.map_id.uuid[0] ^= 0xffu;
            }
            raystar_interfaces::msg::PathResult path_result;
            path_result.cost = 1.0;
            path_result.path.header.frame_id = "map";
            path_result.path.poses.push_back(goal_handle->get_goal()->start);
            path_result.path.poses.push_back(goal_handle->get_goal()->goal);
            result->path_results.push_back(std::move(path_result));
            goal_handle->succeed(result);
            success_count_.fetch_add(1);
          } catch (...) {
            worker_error_count_.fetch_add(1);
          }
        });
      });
    executor_.add_node(node_);
    executor_.add_node(client_node_);
    spin_thread_ = std::thread([this]() { executor_.spin(); });
  }

  ~TestActionServer() {
    {
      std::lock_guard<std::mutex> lock(worker_mutex_);
      stopping_ = true;
    }
    worker_cv_.notify_all();
    executor_.cancel();
    if (spin_thread_.joinable()) {
      spin_thread_.join();
    }
    for (auto& worker : workers_) {
      if (worker.joinable()) {
        worker.join();
      }
    }
    server_.reset();
    node_.reset();
    client_node_.reset();
  }

  rclcpp::Node::SharedPtr clientNode() const {
    return client_node_;
  }
  std::size_t goalCount() const {
    return goal_count_.load();
  }
  std::size_t cancelCount() const {
    return cancel_count_.load();
  }
  std::size_t successCount() const {
    return success_count_.load();
  }
  std::size_t canceledCount() const {
    return canceled_count_.load();
  }
  std::size_t workerErrorCount() const {
    return worker_error_count_.load();
  }
  std::optional<PlanningAction::Goal> lastGoal() const {
    std::lock_guard<std::mutex> lock(goal_mutex_);
    return last_goal_;
  }

  bool waitForClient(std::chrono::milliseconds timeout) {
    auto client = rclcpp_action::create_client<PlanningAction>(client_node_, action_name_);
    return client->wait_for_action_server(timeout);
  }

private:
  static std::string uniqueName(const char* suffix) {
    static std::atomic<unsigned int> sequence{0};
    return std::string("raystar_panel_") + suffix + "_" + std::to_string(sequence.fetch_add(1));
  }

  void finishCanceled(const std::shared_ptr<GoalHandle>& goal_handle) {
    auto result = std::make_shared<PlanningAction::Result>();
    result->success = false;
    result->result_info.status = raystar_interfaces::msg::PlanningResultInfo::STATUS_CANCELLED;
    result->result_info.limits_reached =
      raystar_interfaces::msg::PlanningResultInfo::LIMIT_CANCELLED;
    result->result_info.map_id = goal_handle->get_goal()->map_id;
    goal_handle->canceled(result);
    canceled_count_.fetch_add(1);
  }

  rclcpp::Node::SharedPtr node_;
  rclcpp::Node::SharedPtr client_node_;
  std::string action_name_;
  ServerMode mode_;
  std::chrono::milliseconds execution_delay_;
  rclcpp_action::Server<PlanningAction>::SharedPtr server_;
  rclcpp::executors::MultiThreadedExecutor executor_;
  std::thread spin_thread_;
  mutable std::mutex worker_mutex_;
  std::condition_variable worker_cv_;
  std::vector<std::thread> workers_;
  bool stopping_{false};
  mutable std::mutex goal_mutex_;
  std::optional<PlanningAction::Goal> last_goal_;
  std::atomic<std::size_t> goal_count_{0};
  std::atomic<std::size_t> cancel_count_{0};
  std::atomic<std::size_t> success_count_{0};
  std::atomic<std::size_t> canceled_count_{0};
  std::atomic<std::size_t> worker_error_count_{0};
};

class TestGoalSetActionServer {
public:
  TestGoalSetActionServer(
    std::string action_name,
    ServerMode mode,
    std::chrono::milliseconds execution_delay = std::chrono::milliseconds::zero())
    : node_(rclcpp::Node::make_shared(uniqueName("goal_set_server")))
    , client_node_(rclcpp::Node::make_shared(uniqueName("goal_set_client")))
    , action_name_(std::move(action_name))
    , mode_(mode)
    , execution_delay_(execution_delay) {
    server_ = rclcpp_action::create_server<GoalSetAction>(
      node_,
      action_name_,
      [this](const rclcpp_action::GoalUUID&, std::shared_ptr<const GoalSetAction::Goal> goal) {
        {
          std::lock_guard<std::mutex> lock(goal_mutex_);
          last_goal_ = *goal;
        }
        goal_count_.fetch_add(1);
        return rclcpp_action::GoalResponse::ACCEPT_AND_EXECUTE;
      },
      [this](const std::shared_ptr<GoalSetGoalHandle>&) {
        cancel_request_count_.fetch_add(1);
        return rclcpp_action::CancelResponse::ACCEPT;
      },
      [this](const std::shared_ptr<GoalSetGoalHandle>& goal_handle) {
        std::lock_guard<std::mutex> lock(worker_mutex_);
        if (stopping_) {
          return;
        }
        workers_.emplace_back([this, goal_handle]() noexcept {
          try {
            const auto deadline = std::chrono::steady_clock::now() + execution_delay_;
            while (std::chrono::steady_clock::now() < deadline) {
              if (goal_handle->is_canceling()) {
                finishCanceled(goal_handle);
                return;
              }
              std::unique_lock<std::mutex> lock(worker_mutex_);
              if (worker_cv_.wait_for(
                    lock, std::chrono::milliseconds(5), [this]() { return stopping_; })) {
                return;
              }
            }
            {
              std::lock_guard<std::mutex> lock(worker_mutex_);
              if (stopping_) {
                return;
              }
            }
            if (goal_handle->is_canceling()) {
              finishCanceled(goal_handle);
              return;
            }
            auto result = buildCompleteResult(*goal_handle->get_goal());
            goal_handle->succeed(result);
            terminal_success_count_.fetch_add(1);
          } catch (...) {
            worker_error_count_.fetch_add(1);
          }
        });
      });
    executor_.add_node(node_);
    executor_.add_node(client_node_);
    spin_thread_ = std::thread([this]() { executor_.spin(); });
  }

  ~TestGoalSetActionServer() {
    {
      std::lock_guard<std::mutex> lock(worker_mutex_);
      stopping_ = true;
    }
    worker_cv_.notify_all();
    executor_.cancel();
    if (spin_thread_.joinable()) {
      spin_thread_.join();
    }
    for (auto& worker : workers_) {
      if (worker.joinable()) {
        worker.join();
      }
    }
    server_.reset();
    node_.reset();
    client_node_.reset();
  }

  rclcpp::Node::SharedPtr clientNode() const {
    return client_node_;
  }
  std::size_t goalCount() const {
    return goal_count_.load();
  }
  std::size_t cancelRequestCount() const {
    return cancel_request_count_.load();
  }
  std::size_t terminalSuccessCount() const {
    return terminal_success_count_.load();
  }
  std::size_t terminalCanceledCount() const {
    return terminal_canceled_count_.load();
  }
  std::size_t workerErrorCount() const {
    return worker_error_count_.load();
  }
  std::optional<GoalSetAction::Goal> lastGoal() const {
    std::lock_guard<std::mutex> lock(goal_mutex_);
    return last_goal_;
  }

  bool waitForClient(std::chrono::milliseconds timeout) {
    auto client = rclcpp_action::create_client<GoalSetAction>(client_node_, action_name_);
    return client->wait_for_action_server(timeout);
  }

private:
  static std::string uniqueName(const char* suffix) {
    static std::atomic<unsigned int> sequence{0};
    return std::string("raystar_panel_") + suffix + "_" + std::to_string(sequence.fetch_add(1));
  }

  std::shared_ptr<GoalSetAction::Result> buildCompleteResult(
    const GoalSetAction::Goal& goal) const {
    auto result = std::make_shared<GoalSetAction::Result>();
    result->success = true;
    result->requested_start = goal.start;
    result->requested_allow_self_crossing = goal.allow_self_crossing;
    result->requested_allow_unknown = goal.allow_unknown;

    auto& aggregate = result->result_info;
    aggregate.map_id = goal.map_id;
    aggregate.environment_id.uuid[0] = 0x42u;
    aggregate.status = raystar_interfaces::msg::PlanningResultInfo::STATUS_COMPLETE;
    aggregate.limits_reached = raystar_interfaces::msg::PlanningResultInfo::LIMIT_NONE;
    aggregate.request_satisfied = true;
    aggregate.search_complete = true;
    aggregate.output_complete = true;
    aggregate.debug_requested = goal.include_debug;
    aggregate.debug_output_complete = true;
    aggregate.requested_goal_count = static_cast<std::uint32_t>(goal.goals.size());
    aggregate.returned_goal_count = static_cast<std::uint32_t>(goal.goals.size());
    aggregate.completed_goal_count = static_cast<std::uint32_t>(goal.goals.size());

    result->goal_results.reserve(goal.goals.size());
    for (std::size_t index = 0; index < goal.goals.size(); ++index) {
      raystar_interfaces::msg::GoalPathResult goal_result;
      goal_result.goal_index = static_cast<std::uint32_t>(index);
      goal_result.goal = goal.goals[index];
      goal_result.requested_max_path_length = goal.max_path_lengths[index];
      auto& info = goal_result.result_info;
      info.map_id = goal.map_id;
      info.environment_id = aggregate.environment_id;
      info.limits_reached = raystar_interfaces::msg::PlanningResultInfo::LIMIT_NONE;
      info.request_satisfied = true;
      info.search_complete = true;
      info.output_complete = true;
      info.debug_requested = goal.include_debug;
      info.debug_output_complete = true;
      info.search_mode = PlanningAction::Goal::SEARCH_MODE_ALL_WITHIN_LENGTH;
      info.requested_max_path_length = goal.max_path_lengths[index];
      info.cost_bound_exhausted = true;

      // Keep the middle goal as a certified empty result.  This deliberately
      // exercises the distinction between payload presence and bounded-search
      // completeness, including when duplicate target coordinates are used.
      const bool has_path = index != 1u;
      goal_result.success = has_path;
      info.status = has_path ? raystar_interfaces::msg::PlanningResultInfo::STATUS_COMPLETE
                             : raystar_interfaces::msg::PlanningResultInfo::STATUS_NO_PATH;
      if (has_path) {
        raystar_interfaces::msg::PathResult path_result;
        path_result.cost = std::min(goal.max_path_lengths[index], static_cast<double>(index + 1u));
        path_result.path.header.frame_id = goal.goals[index].header.frame_id;
        path_result.path.poses.push_back(goal.start);
        path_result.path.poses.push_back(goal.goals[index]);
        goal_result.path_results.push_back(std::move(path_result));
        info.found_path_count = 1;
        info.returned_path_count = 1;
        ++aggregate.goals_with_paths;
        ++aggregate.found_path_count;
        ++aggregate.returned_path_count;
      } else {
        goal_result.message = "No path exists within this goal's inclusive bound.";
      }
      result->goal_results.emplace_back(std::move(goal_result));
    }

    if (mode_ == ServerMode::mismatched_map) {
      aggregate.map_id.uuid[0] ^= 0xffu;
    } else if (mode_ == ServerMode::inner_bound_mismatch && !result->goal_results.empty()) {
      // The outer association remains correct.  Only the structured
      // PlanningResultInfo echo is corrupt, so a client that validates just
      // GoalPathResult.requested_max_path_length would accept stale metadata.
      result->goal_results.front().result_info.requested_max_path_length += 1.0;
    } else if (mode_ == ServerMode::partial_output && !result->goal_results.empty()) {
      // A partial response is still a structurally valid Action result and
      // may contain useful paths.  response.success describes aggregate
      // request satisfaction; per-goal success continues to describe payload
      // presence only.
      result->success = false;
      aggregate.status = raystar_interfaces::msg::PlanningResultInfo::STATUS_PARTIAL_OUTPUT;
      aggregate.limits_reached = raystar_interfaces::msg::PlanningResultInfo::LIMIT_MAX_PATHS;
      aggregate.request_satisfied = false;
      aggregate.output_complete = false;
      aggregate.completed_goal_count =
        aggregate.requested_goal_count > 0u ? aggregate.requested_goal_count - 1u : 0u;
      ++aggregate.found_path_count;

      auto& first = result->goal_results.front();
      first.result_info.status = raystar_interfaces::msg::PlanningResultInfo::STATUS_PARTIAL_OUTPUT;
      first.result_info.limits_reached =
        raystar_interfaces::msg::PlanningResultInfo::LIMIT_MAX_PATHS;
      first.result_info.request_satisfied = false;
      first.result_info.output_complete = false;
      first.result_info.found_path_count = 2;
    }
    result->message = "Synthetic mixed complete/empty goal-set result.";
    return result;
  }

  void finishCanceled(const std::shared_ptr<GoalSetGoalHandle>& goal_handle) {
    auto result = std::make_shared<GoalSetAction::Result>();
    result->success = false;
    result->requested_start = goal_handle->get_goal()->start;
    result->requested_allow_self_crossing = goal_handle->get_goal()->allow_self_crossing;
    result->requested_allow_unknown = goal_handle->get_goal()->allow_unknown;
    result->result_info.map_id = goal_handle->get_goal()->map_id;
    result->result_info.status = raystar_interfaces::msg::PlanningResultInfo::STATUS_CANCELLED;
    result->result_info.limits_reached =
      raystar_interfaces::msg::PlanningResultInfo::LIMIT_CANCELLED;
    result->result_info.requested_goal_count =
      static_cast<std::uint32_t>(goal_handle->get_goal()->goals.size());
    goal_handle->canceled(result);
    terminal_canceled_count_.fetch_add(1);
  }

  rclcpp::Node::SharedPtr node_;
  rclcpp::Node::SharedPtr client_node_;
  std::string action_name_;
  ServerMode mode_;
  std::chrono::milliseconds execution_delay_;
  rclcpp_action::Server<GoalSetAction>::SharedPtr server_;
  rclcpp::executors::MultiThreadedExecutor executor_;
  std::thread spin_thread_;
  mutable std::mutex worker_mutex_;
  std::condition_variable worker_cv_;
  std::vector<std::thread> workers_;
  bool stopping_{false};
  mutable std::mutex goal_mutex_;
  std::optional<GoalSetAction::Goal> last_goal_;
  std::atomic<std::size_t> goal_count_{0};
  std::atomic<std::size_t> cancel_request_count_{0};
  std::atomic<std::size_t> terminal_success_count_{0};
  std::atomic<std::size_t> terminal_canceled_count_{0};
  std::atomic<std::size_t> worker_error_count_{0};
};

class TestMapPublisher {
public:
  explicit TestMapPublisher(rclcpp::Node::SharedPtr node, std::string topic)
    : node_(std::move(node)), topic_(std::move(topic)) {
    publisher_ = node_->create_publisher<nav_msgs::msg::OccupancyGrid>(
      topic_, rclcpp::QoS(1).transient_local().reliable());
  }

  nav_msgs::msg::OccupancyGrid publishMap() {
    nav_msgs::msg::OccupancyGrid map;
    map.header.frame_id = "map";
    map.info.width = 4;
    map.info.height = 4;
    map.info.resolution = 1.0f;
    map.info.origin.orientation.w = 1.0;
    map.data.assign(map.info.width * map.info.height, 0);
    publisher_->publish(map);
    return map;
  }

  void publish(const nav_msgs::msg::OccupancyGrid& map) {
    publisher_->publish(map);
  }

private:
  rclcpp::Node::SharedPtr node_;
  std::string topic_;
  rclcpp::Publisher<nav_msgs::msg::OccupancyGrid>::SharedPtr publisher_;
};

class TestClickedPointPublisher {
public:
  TestClickedPointPublisher(rclcpp::Node::SharedPtr node, std::string topic)
    : node_(std::move(node)), topic_(std::move(topic)) {
    publisher_ =
      node_->create_publisher<geometry_msgs::msg::PointStamped>(topic_, rclcpp::QoS(10).reliable());
  }

  void publish(double x, double y, const std::string& frame_id = "map") {
    geometry_msgs::msg::PointStamped point;
    point.header.frame_id = frame_id;
    point.point.x = x;
    point.point.y = y;
    publisher_->publish(point);
  }

  std::size_t subscriptionCount() const {
    return publisher_->get_subscription_count();
  }

private:
  rclcpp::Node::SharedPtr node_;
  std::string topic_;
  rclcpp::Publisher<geometry_msgs::msg::PointStamped>::SharedPtr publisher_;
};

template <typename Predicate>
bool waitForGui(Predicate predicate, std::chrono::milliseconds timeout) {
  const auto deadline = std::chrono::steady_clock::now() + timeout;
  while (std::chrono::steady_clock::now() < deadline) {
    QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
    if (predicate()) {
      return true;
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
  }
  QCoreApplication::processEvents(QEventLoop::AllEvents, 10);
  return predicate();
}

std::string nextTestName(const char* stem) {
  static std::atomic<unsigned int> sequence{0};
  return std::string("/") + stem + "_" + std::to_string(sequence.fetch_add(1));
}

std::unique_ptr<raystar_rviz_plugins::RaystarPanel> makeInitializedPanel(
  TestDisplayContext& context,
  const std::string& action_name,
  const std::string& map_topic,
  std::chrono::milliseconds timeout,
  const std::string& goal_set_action_name = {},
  const std::string& clicked_point_topic = {}) {
  auto panel = std::make_unique<raystar_rviz_plugins::RaystarPanel>(nullptr, timeout);
  panel->findChild<QLineEdit*>("action_name_edit")->setText(QString::fromStdString(action_name));
  panel->findChild<QLineEdit*>("map_topic_edit")->setText(QString::fromStdString(map_topic));
  if (!goal_set_action_name.empty()) {
    panel->findChild<QLineEdit*>("goal_set_action_name_edit")
      ->setText(QString::fromStdString(goal_set_action_name));
  }
  if (!clicked_point_topic.empty()) {
    panel->findChild<QLineEdit*>("clicked_point_topic_edit")
      ->setText(QString::fromStdString(clicked_point_topic));
  }
  panel->initialize(&context);
  return panel;
}

QLabel* statusLabel(raystar_rviz_plugins::RaystarPanel& panel) {
  return panel.findChild<QLabel*>("status_label");
}

QPushButton* planButton(raystar_rviz_plugins::RaystarPanel& panel) {
  return panel.findChild<QPushButton*>("plan_button");
}

QTableWidget* resultsTable(raystar_rviz_plugins::RaystarPanel& panel) {
  return panel.findChild<QTableWidget*>("results_table");
}

QTableWidget* goalResultsTable(raystar_rviz_plugins::RaystarPanel& panel) {
  return panel.findChild<QTableWidget*>("goal_results_table");
}

QTableWidget* multiGoalTable(raystar_rviz_plugins::RaystarPanel& panel) {
  return panel.findChild<QTableWidget*>("multi_goal_table");
}

QComboBox* searchModeCombo(raystar_rviz_plugins::RaystarPanel& panel) {
  return panel.findChild<QComboBox*>("search_mode_combo");
}

QDoubleSpinBox* maxPathLengthSpinbox(raystar_rviz_plugins::RaystarPanel& panel) {
  return panel.findChild<QDoubleSpinBox*>("max_path_length_spinbox");
}

struct TestGoalRow {
  double x;
  double y;
  double max_path_length;
};

void setMultiGoalRows(raystar_rviz_plugins::RaystarPanel& panel,
                      const std::vector<TestGoalRow>& goals) {
  auto* table = multiGoalTable(panel);
  ASSERT_NE(nullptr, table);
  table->setRowCount(static_cast<int>(goals.size()));
  for (std::size_t index = 0; index < goals.size(); ++index) {
    const int row = static_cast<int>(index);
    table->setItem(row, 0, new QTableWidgetItem(QString::number(goals[index].x, 'g', 17)));
    table->setItem(row, 1, new QTableWidgetItem(QString::number(goals[index].y, 'g', 17)));
    table->setItem(
      row, 2, new QTableWidgetItem(QString::number(goals[index].max_path_length, 'g', 17)));
  }
}

QString tableText(const QTableWidget& table, int row, int column) {
  const auto* item = table.item(row, column);
  return item ? item->text() : QString();
}

bool waitForMapAndPlanReady(raystar_rviz_plugins::RaystarPanel& panel,
                            TestMapPublisher& publisher,
                            const nav_msgs::msg::OccupancyGrid& map,
                            std::chrono::milliseconds timeout) {
  auto* label = statusLabel(panel);
  auto* button = planButton(panel);
  if (!label || !button) {
    return false;
  }
  return waitForGui(
    [&]() {
      if (!label->text().contains(QStringLiteral("Map received"))) {
        publisher.publish(map);
      }
      return label->text().contains(QStringLiteral("Map received")) && button->isEnabled();
    },
    timeout);
}

constexpr char kPluginId[] = "raystar_rviz_plugins/RaystarPanel";
constexpr char kLegacyPluginId[] = "raystar_rviz_plugins::RaystarPanel";
constexpr char kPluginType[] = "raystar_rviz_plugins::RaystarPanel";

using PanelLoader = pluginlib::ClassLoader<rviz_common::Panel>;

std::shared_ptr<rviz_common::Panel> loadPanel(PanelLoader& loader) {
  return loader.createSharedInstance(kPluginId);
}

TEST(RaystarPanelPlugin, DeclaresOneStableIdAndLoads) {
  PanelLoader loader("rviz_common", "rviz_common::Panel");
  const auto declared = loader.getDeclaredClasses();

  EXPECT_EQ(1, std::count(declared.begin(), declared.end(), kPluginId));
  EXPECT_EQ(0, std::count(declared.begin(), declared.end(), kLegacyPluginId));
  ASSERT_TRUE(loader.isClassAvailable(kPluginId));
  EXPECT_EQ(kPluginType, loader.getClassType(kPluginId));

  auto panel = loadPanel(loader);
  ASSERT_NE(nullptr, panel);
  panel.reset();
}

TEST(RaystarPanelPlugin, PersistsConnectionAndPlanningOptions) {
  PanelLoader loader("rviz_common", "rviz_common::Panel");
  auto panel = loadPanel(loader);
  ASSERT_NE(nullptr, panel);

  rviz_common::Config input;
  input.mapSetValue("action_name", QString("/robot1/raystar/plan_paths"));
  input.mapSetValue("map_topic", QString("/robot1/map"));
  input.mapSetValue("start_x", QString("1.25"));
  input.mapSetValue("start_y", QString("2.5"));
  input.mapSetValue("goal_x", QString("3.75"));
  input.mapSetValue("goal_y", QString("4.5"));
  input.mapSetValue("search_mode", 1);
  input.mapSetValue("k", 7);
  input.mapSetValue("max_path_length", QString("42.125"));
  input.mapSetValue("allow_self_crossing", true);
  input.mapSetValue("allow_unknown", true);
  input.mapSetValue("request_debug", true);
  panel->load(input);

  rviz_common::Config output;
  panel->save(output);

  QString text;
  int integer = 0;
  bool boolean = false;
  ASSERT_TRUE(output.mapGetString("action_name", &text));
  EXPECT_EQ(QString("/robot1/raystar/plan_paths"), text);
  ASSERT_TRUE(output.mapGetString("map_topic", &text));
  EXPECT_EQ(QString("/robot1/map"), text);
  ASSERT_TRUE(output.mapGetString("start_x", &text));
  EXPECT_EQ(QString("1.25"), text);
  ASSERT_TRUE(output.mapGetString("start_y", &text));
  EXPECT_EQ(QString("2.5"), text);
  ASSERT_TRUE(output.mapGetString("goal_x", &text));
  EXPECT_EQ(QString("3.75"), text);
  ASSERT_TRUE(output.mapGetString("goal_y", &text));
  EXPECT_EQ(QString("4.5"), text);
  ASSERT_TRUE(output.mapGetInt("k", &integer));
  EXPECT_EQ(7, integer);
  ASSERT_TRUE(output.mapGetInt("search_mode", &integer));
  EXPECT_EQ(1, integer);
  ASSERT_TRUE(output.mapGetString("max_path_length", &text));
  EXPECT_DOUBLE_EQ(42.125, text.toDouble());
  ASSERT_TRUE(output.mapGetBool("allow_self_crossing", &boolean));
  EXPECT_TRUE(boolean);
  ASSERT_TRUE(output.mapGetBool("allow_unknown", &boolean));
  EXPECT_TRUE(boolean);
  ASSERT_TRUE(output.mapGetBool("request_debug", &boolean));
  EXPECT_TRUE(boolean);

  panel.reset();
}

TEST(RaystarPanelPlugin, PersistsMultiGoalListOfMapsAndEndpoints) {
  PanelLoader loader("rviz_common", "rviz_common::Panel");
  auto panel = loadPanel(loader);
  ASSERT_NE(nullptr, panel);

  const std::vector<TestGoalRow> expected_goals = {
    {3.25, 4.5, 17.5},
    {3.25, 4.5, 23.75},
    {-2.0, 6.125, 31.125},
  };
  rviz_common::Config input;
  input.mapSetValue("goal_set_action_name", QString("/robot1/raystar/plan_goal_set"));
  input.mapSetValue("clicked_point_topic", QString("/robot1/clicked_point"));
  input.mapSetValue("search_mode", 2);
  auto input_goals = input.mapMakeChild("multi_goals");
  for (const auto& expected : expected_goals) {
    auto item = input_goals.listAppendNew();
    item.mapSetValue("x", QString::number(expected.x, 'g', 17));
    item.mapSetValue("y", QString::number(expected.y, 'g', 17));
    item.mapSetValue("max_path_length", QString::number(expected.max_path_length, 'g', 17));
  }
  panel->load(input);

  auto* editor = panel->findChild<QTableWidget*>("multi_goal_table");
  ASSERT_NE(nullptr, editor);
  ASSERT_EQ(static_cast<int>(expected_goals.size()), editor->rowCount());
  for (std::size_t index = 0; index < expected_goals.size(); ++index) {
    const int row = static_cast<int>(index);
    EXPECT_DOUBLE_EQ(expected_goals[index].x, tableText(*editor, row, 0).toDouble());
    EXPECT_DOUBLE_EQ(expected_goals[index].y, tableText(*editor, row, 1).toDouble());
    EXPECT_DOUBLE_EQ(expected_goals[index].max_path_length, tableText(*editor, row, 2).toDouble());
  }

  rviz_common::Config output;
  panel->save(output);

  QString text;
  int search_mode = -1;
  ASSERT_TRUE(output.mapGetString("goal_set_action_name", &text));
  EXPECT_EQ(QString("/robot1/raystar/plan_goal_set"), text);
  ASSERT_TRUE(output.mapGetString("clicked_point_topic", &text));
  EXPECT_EQ(QString("/robot1/clicked_point"), text);
  ASSERT_TRUE(output.mapGetInt("search_mode", &search_mode));
  EXPECT_EQ(2, search_mode);
  const auto output_goals = output.mapGetChild("multi_goals");
  ASSERT_TRUE(output_goals.isValid());
  ASSERT_EQ(static_cast<int>(expected_goals.size()), output_goals.listLength());
  for (std::size_t index = 0; index < expected_goals.size(); ++index) {
    const auto item = output_goals.listChildAt(static_cast<int>(index));
    ASSERT_TRUE(item.mapGetString("x", &text));
    EXPECT_DOUBLE_EQ(expected_goals[index].x, text.toDouble());
    ASSERT_TRUE(item.mapGetString("y", &text));
    EXPECT_DOUBLE_EQ(expected_goals[index].y, text.toDouble());
    ASSERT_TRUE(item.mapGetString("max_path_length", &text));
    EXPECT_DOUBLE_EQ(expected_goals[index].max_path_length, text.toDouble());
  }

  panel.reset();
}

TEST(RaystarPanelPlugin, OlderConfigsKeepSafeDefaultsForNewFields) {
  PanelLoader loader("rviz_common", "rviz_common::Panel");
  auto panel = loadPanel(loader);
  ASSERT_NE(nullptr, panel);

  rviz_common::Config old_config;
  old_config.mapSetValue("map_topic", QString("/legacy_map"));
  panel->load(old_config);

  rviz_common::Config output;
  panel->save(output);

  QString action_name;
  QString goal_set_action_name;
  QString clicked_point_topic;
  QString max_path_length;
  int search_mode = -1;
  bool boolean = true;
  ASSERT_TRUE(output.mapGetString("action_name", &action_name));
  EXPECT_EQ(QString("/raystar/plan_paths"), action_name);
  ASSERT_TRUE(output.mapGetString("goal_set_action_name", &goal_set_action_name));
  EXPECT_EQ(QString("/raystar/plan_goal_set"), goal_set_action_name);
  ASSERT_TRUE(output.mapGetString("clicked_point_topic", &clicked_point_topic));
  EXPECT_EQ(QString("/clicked_point"), clicked_point_topic);
  ASSERT_TRUE(output.mapGetInt("search_mode", &search_mode));
  EXPECT_EQ(0, search_mode);
  ASSERT_TRUE(output.mapGetString("max_path_length", &max_path_length));
  EXPECT_GT(max_path_length.toDouble(), 0.0);
  ASSERT_TRUE(output.mapGetBool("allow_self_crossing", &boolean));
  EXPECT_FALSE(boolean);
  ASSERT_TRUE(output.mapGetBool("allow_unknown", &boolean));
  EXPECT_FALSE(boolean);
  ASSERT_TRUE(output.mapGetBool("request_debug", &boolean));
  EXPECT_FALSE(boolean);
  const auto goals = output.mapGetChild("multi_goals");
  ASSERT_TRUE(goals.isValid());
  ASSERT_EQ(1, goals.listLength());

  panel.reset();
}

TEST(RaystarPanelInteraction, InitializesWithDisplayContextAndRendersActionResult) {
  const auto action_name = nextTestName("plan_paths");
  const auto map_topic = nextTestName("map");
  TestActionServer server(action_name, ServerMode::success);
  ASSERT_TRUE(server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(server.clientNode(), map_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(context, action_name, map_topic, std::chrono::seconds(2));

  auto map = map_publisher.publishMap();
  auto* label = statusLabel(*panel);
  auto* button = planButton(*panel);
  auto* table = resultsTable(*panel);
  ASSERT_NE(nullptr, label);
  ASSERT_NE(nullptr, button);
  ASSERT_NE(nullptr, table);
  ASSERT_TRUE(waitForGui(
    [&]() {
      if (!label->text().contains(QStringLiteral("Map received"))) {
        map_publisher.publish(map);
      }
      return label->text().contains(QStringLiteral("Map received"));
    },
    std::chrono::seconds(2)))
    << label->text().toStdString();
  ASSERT_TRUE(waitForGui([&]() { return button->isEnabled(); }, std::chrono::seconds(2)))
    << label->text().toStdString();

  button->click();
  ASSERT_TRUE(waitForGui([&]() { return label->text().contains(QStringLiteral("SUCCEEDED")); },
                         std::chrono::seconds(2)))
    << label->text().toStdString();
  EXPECT_EQ(1, table->rowCount());
  EXPECT_EQ(1u, server.goalCount());
  EXPECT_EQ(1u, server.successCount());
}

TEST(RaystarPanelInteraction, RejectsResultForDifferentMapId) {
  const auto action_name = nextTestName("plan_paths");
  const auto map_topic = nextTestName("map");
  TestActionServer server(action_name, ServerMode::mismatched_map);
  ASSERT_TRUE(server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(server.clientNode(), map_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(context, action_name, map_topic, std::chrono::seconds(2));
  auto map = map_publisher.publishMap();
  auto* label = statusLabel(*panel);
  auto* button = planButton(*panel);
  ASSERT_TRUE(waitForGui(
    [&]() {
      if (!label->text().contains(QStringLiteral("Map received"))) {
        map_publisher.publish(map);
      }
      return label->text().contains(QStringLiteral("Map received"));
    },
    std::chrono::seconds(2)));
  ASSERT_TRUE(waitForGui([&]() { return button->isEnabled(); }, std::chrono::seconds(2)));

  button->click();
  ASSERT_TRUE(
    waitForGui([&]() { return label->text().contains(QStringLiteral("different cached map")); },
               std::chrono::seconds(2)))
    << label->text().toStdString();
}

TEST(RaystarPanelInteraction, SendsCostBoundedEnumerationGoal) {
  const auto action_name = nextTestName("plan_paths");
  const auto map_topic = nextTestName("map");
  TestActionServer server(action_name, ServerMode::success);
  ASSERT_TRUE(server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(server.clientNode(), map_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(context, action_name, map_topic, std::chrono::seconds(2));

  auto* mode = searchModeCombo(*panel);
  auto* bound = maxPathLengthSpinbox(*panel);
  auto* label = statusLabel(*panel);
  auto* button = planButton(*panel);
  ASSERT_NE(nullptr, mode);
  ASSERT_NE(nullptr, bound);
  ASSERT_NE(nullptr, label);
  ASSERT_NE(nullptr, button);
  mode->setCurrentIndex(1);
  bound->setValue(37.25);
  EXPECT_FALSE(panel->findChild<QSpinBox*>("k_spinbox")->isEnabled());
  EXPECT_TRUE(bound->isEnabled());

  auto map = map_publisher.publishMap();
  ASSERT_TRUE(waitForGui(
    [&]() {
      if (!label->text().contains(QStringLiteral("Map received"))) {
        map_publisher.publish(map);
      }
      return label->text().contains(QStringLiteral("Map received"));
    },
    std::chrono::seconds(2)));
  ASSERT_TRUE(waitForGui([&]() { return button->isEnabled(); }, std::chrono::seconds(2)));

  button->click();
  ASSERT_TRUE(waitForGui([&]() { return server.lastGoal().has_value(); }, std::chrono::seconds(2)));
  const auto goal = server.lastGoal();
  ASSERT_TRUE(goal.has_value());
  EXPECT_EQ(goal->search_mode, PlanningAction::Goal::SEARCH_MODE_ALL_WITHIN_LENGTH);
  EXPECT_EQ(goal->k, 0);
  EXPECT_DOUBLE_EQ(goal->max_path_length, 37.25);
  ASSERT_TRUE(waitForGui([&]() { return server.successCount() == 1u; }, std::chrono::seconds(2)));
  EXPECT_EQ(0u, server.workerErrorCount());
}

TEST(RaystarPanelInteraction, CancelsShortTimeoutGoal) {
  const auto action_name = nextTestName("plan_paths");
  const auto map_topic = nextTestName("map");
  TestActionServer server(action_name, ServerMode::delayed, std::chrono::seconds(2));
  ASSERT_TRUE(server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(server.clientNode(), map_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(context, action_name, map_topic, std::chrono::milliseconds(80));
  auto map = map_publisher.publishMap();
  auto* label = statusLabel(*panel);
  auto* button = planButton(*panel);
  ASSERT_TRUE(waitForGui(
    [&]() {
      if (!label->text().contains(QStringLiteral("Map received"))) {
        map_publisher.publish(map);
      }
      return label->text().contains(QStringLiteral("Map received"));
    },
    std::chrono::seconds(2)));
  ASSERT_TRUE(waitForGui([&]() { return button->isEnabled(); }, std::chrono::seconds(2)));

  button->click();
  ASSERT_TRUE(waitForGui(
    [&]() { return label->text().contains(QStringLiteral("timed out and was canceled")); },
    std::chrono::seconds(2)))
    << label->text().toStdString();
  ASSERT_TRUE(waitForGui([&]() { return server.cancelCount() > 0; }, std::chrono::seconds(2)));
  ASSERT_TRUE(waitForGui([&]() { return server.canceledCount() > 0; }, std::chrono::seconds(2)));
  EXPECT_EQ(0u, server.workerErrorCount());
}

TEST(RaystarPanelInteraction, DestructionCancelsGoalAndLeavesNoQueuedUiCallback) {
  const auto action_name = nextTestName("plan_paths");
  const auto map_topic = nextTestName("map");
  TestActionServer server(action_name, ServerMode::delayed, std::chrono::seconds(2));
  ASSERT_TRUE(server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(server.clientNode(), map_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(context, action_name, map_topic, std::chrono::seconds(5));
  auto map = map_publisher.publishMap();
  auto* label = statusLabel(*panel);
  auto* button = planButton(*panel);
  ASSERT_TRUE(waitForGui(
    [&]() {
      if (!label->text().contains(QStringLiteral("Map received"))) {
        map_publisher.publish(map);
      }
      return label->text().contains(QStringLiteral("Map received"));
    },
    std::chrono::seconds(2)));
  ASSERT_TRUE(waitForGui([&]() { return button->isEnabled(); }, std::chrono::seconds(2)));

  button->click();
  ASSERT_TRUE(waitForGui([&]() { return server.goalCount() > 0; }, std::chrono::seconds(2)));
  panel.reset();

  // The action acceptance/result callbacks may run after the QWidget has
  // gone away. The panel must retain only non-UI state long enough to cancel
  // the accepted goal, and must not dereference the destroyed object.
  ASSERT_TRUE(waitForGui([&]() { return server.cancelCount() > 0; }, std::chrono::seconds(2)));
  ASSERT_TRUE(waitForGui([&]() { return server.canceledCount() > 0; }, std::chrono::seconds(2)));
  EXPECT_EQ(0u, server.workerErrorCount());
}

TEST(RaystarPanelInteraction, SendsOneSharedTreeGoalSetWithIndependentBoundsAndRendersMixedResult) {
  const auto action_name = nextTestName("plan_paths");
  const auto goal_set_action_name = nextTestName("plan_goal_set");
  const auto map_topic = nextTestName("map");
  TestActionServer single_server(action_name, ServerMode::success);
  TestGoalSetActionServer goal_set_server(goal_set_action_name, ServerMode::success);
  ASSERT_TRUE(single_server.waitForClient(std::chrono::seconds(2)));
  ASSERT_TRUE(goal_set_server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(goal_set_server.clientNode(), map_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(goal_set_server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(
    context, action_name, map_topic, std::chrono::seconds(2), goal_set_action_name);

  auto* mode = searchModeCombo(*panel);
  auto* label = statusLabel(*panel);
  auto* button = planButton(*panel);
  auto* path_table = resultsTable(*panel);
  auto* goal_table = goalResultsTable(*panel);
  ASSERT_NE(nullptr, mode);
  ASSERT_NE(nullptr, label);
  ASSERT_NE(nullptr, button);
  ASSERT_NE(nullptr, path_table);
  ASSERT_NE(nullptr, goal_table);
  mode->setCurrentIndex(2);
  panel->findChild<QLineEdit*>("start_x_edit")->setText("1.25");
  panel->findChild<QLineEdit*>("start_y_edit")->setText("2.5");
  panel->findChild<QCheckBox*>("allow_self_crossing_checkbox")->setChecked(true);
  panel->findChild<QCheckBox*>("allow_unknown_checkbox")->setChecked(true);
  panel->findChild<QCheckBox*>("request_debug_checkbox")->setChecked(true);

  const std::vector<TestGoalRow> expected_goals = {
    {3.25, 4.5, 17.5},
    // A duplicated target with a different bound must remain a distinct
    // request entry, associated by its input index rather than coordinates.
    {3.25, 4.5, 23.75},
    {-2.0, 6.125, 31.125},
  };
  setMultiGoalRows(*panel, expected_goals);

  const auto map = map_publisher.publishMap();
  ASSERT_TRUE(waitForMapAndPlanReady(*panel, map_publisher, map, std::chrono::seconds(2)))
    << label->text().toStdString();

  button->click();
  ASSERT_TRUE(
    waitForGui([&]() { return goal_set_server.lastGoal().has_value(); }, std::chrono::seconds(2)));
  const auto goal = goal_set_server.lastGoal();
  ASSERT_TRUE(goal.has_value());
  EXPECT_EQ(1u, goal_set_server.goalCount());
  EXPECT_EQ(0u, single_server.goalCount());
  EXPECT_TRUE(raystar_interfaces::mapIdsEqual(goal->map_id, raystar_interfaces::computeMapId(map)));
  EXPECT_EQ("map", goal->start.header.frame_id);
  EXPECT_DOUBLE_EQ(1.25, goal->start.pose.position.x);
  EXPECT_DOUBLE_EQ(2.5, goal->start.pose.position.y);
  EXPECT_DOUBLE_EQ(1.0, goal->start.pose.orientation.w);
  ASSERT_EQ(expected_goals.size(), goal->goals.size());
  ASSERT_EQ(expected_goals.size(), goal->max_path_lengths.size());
  for (std::size_t index = 0; index < expected_goals.size(); ++index) {
    EXPECT_EQ("map", goal->goals[index].header.frame_id);
    EXPECT_DOUBLE_EQ(expected_goals[index].x, goal->goals[index].pose.position.x);
    EXPECT_DOUBLE_EQ(expected_goals[index].y, goal->goals[index].pose.position.y);
    EXPECT_DOUBLE_EQ(1.0, goal->goals[index].pose.orientation.w);
    EXPECT_DOUBLE_EQ(expected_goals[index].max_path_length, goal->max_path_lengths[index]);
  }
  EXPECT_EQ(goal->goals[0], goal->goals[1]);
  EXPECT_NE(goal->max_path_lengths[0], goal->max_path_lengths[1]);
  EXPECT_TRUE(goal->allow_self_crossing);
  EXPECT_TRUE(goal->allow_unknown);
  EXPECT_TRUE(goal->include_debug);

  ASSERT_TRUE(
    waitForGui([&]() { return label->text().contains("SUCCEEDED"); }, std::chrono::seconds(2)))
    << label->text().toStdString();
  ASSERT_EQ(3, goal_table->rowCount());
  EXPECT_EQ(QString("0"), tableText(*goal_table, 0, 0));
  EXPECT_EQ(QString("Complete"), tableText(*goal_table, 0, 4));
  EXPECT_EQ(QString("yes"), tableText(*goal_table, 0, 5));
  EXPECT_EQ(QString("1 / 1"), tableText(*goal_table, 0, 6));
  EXPECT_EQ(QString("1"), tableText(*goal_table, 1, 0));
  EXPECT_EQ(QString("No path"), tableText(*goal_table, 1, 4));
  EXPECT_EQ(QString("yes"), tableText(*goal_table, 1, 5));
  EXPECT_EQ(QString("0 / 0"), tableText(*goal_table, 1, 6));
  EXPECT_EQ(QString("2"), tableText(*goal_table, 2, 0));
  EXPECT_EQ(QString("Complete"), tableText(*goal_table, 2, 4));
  EXPECT_EQ(QString("yes"), tableText(*goal_table, 2, 5));
  EXPECT_EQ(QString("1 / 1"), tableText(*goal_table, 2, 6));

  ASSERT_EQ(2, path_table->rowCount());
  EXPECT_EQ(QString("0"), tableText(*path_table, 0, 0));
  EXPECT_EQ(QString("1"), tableText(*path_table, 0, 1));
  EXPECT_EQ(QString("path_1"), tableText(*path_table, 0, 2));
  EXPECT_DOUBLE_EQ(1.0, tableText(*path_table, 0, 3).toDouble());
  EXPECT_EQ(QString("2"), tableText(*path_table, 1, 0));
  EXPECT_EQ(QString("1"), tableText(*path_table, 1, 1));
  EXPECT_EQ(QString("path_2"), tableText(*path_table, 1, 2));
  EXPECT_DOUBLE_EQ(3.0, tableText(*path_table, 1, 3).toDouble());
  EXPECT_TRUE(label->text().contains("completed 3 of 3 goals"));
  EXPECT_TRUE(label->text().contains("2 goal(s) have paths"));
  EXPECT_TRUE(label->text().contains("Request satisfied: yes"));
  EXPECT_TRUE(label->text().contains("path output complete: yes"));
  EXPECT_FALSE(label->text().contains("Warning:"));
  ASSERT_TRUE(waitForGui([&]() { return goal_set_server.terminalSuccessCount() == 1u; },
                         std::chrono::seconds(2)));
  EXPECT_EQ(0u, goal_set_server.workerErrorCount());
}

TEST(RaystarPanelInteraction, ChangingMultiGoalBulkBudgetUpdatesEveryExistingGoalAndEveryRequest) {
  const auto action_name = nextTestName("plan_paths_unused");
  const auto goal_set_action_name = nextTestName("plan_goal_set");
  const auto map_topic = nextTestName("map");
  TestActionServer single_server(action_name, ServerMode::success);
  TestGoalSetActionServer goal_set_server(goal_set_action_name, ServerMode::success);
  ASSERT_TRUE(single_server.waitForClient(std::chrono::seconds(2)));
  ASSERT_TRUE(goal_set_server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(goal_set_server.clientNode(), map_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(goal_set_server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(
    context, action_name, map_topic, std::chrono::seconds(2), goal_set_action_name);

  auto* mode = searchModeCombo(*panel);
  auto* bulk_budget = maxPathLengthSpinbox(*panel);
  auto* bulk_budget_label = panel->findChild<QLabel*>("max_path_length_label");
  auto* input_table = multiGoalTable(*panel);
  auto* result_table = resultsTable(*panel);
  auto* goal_result_table = goalResultsTable(*panel);
  auto* label = statusLabel(*panel);
  auto* button = planButton(*panel);
  ASSERT_NE(nullptr, mode);
  ASSERT_NE(nullptr, bulk_budget);
  ASSERT_NE(nullptr, bulk_budget_label);
  ASSERT_NE(nullptr, input_table);
  ASSERT_NE(nullptr, result_table);
  ASSERT_NE(nullptr, goal_result_table);
  ASSERT_NE(nullptr, label);
  ASSERT_NE(nullptr, button);

  mode->setCurrentIndex(2);
  EXPECT_EQ(QString("Set all budgets (m):"), bulk_budget_label->text());
  setMultiGoalRows(*panel, {{3.25, 4.5, 7.0}, {-2.0, 6.125, 6.0}});

  const auto map = map_publisher.publishMap();
  ASSERT_TRUE(waitForMapAndPlanReady(*panel, map_publisher, map, std::chrono::seconds(2)))
    << label->text().toStdString();

  const std::vector<double> decreasing_bounds = {7.0, 6.0, 5.0};
  std::size_t expected_goal_count = 0;
  for (const double bound : decreasing_bounds) {
    bulk_budget->setValue(bound);
    for (int row = 0; row < input_table->rowCount(); ++row) {
      EXPECT_DOUBLE_EQ(bound, tableText(*input_table, row, 2).toDouble());
    }
    EXPECT_EQ(0, result_table->rowCount());
    EXPECT_EQ(0, goal_result_table->rowCount());
    EXPECT_TRUE(label->text().contains("Press Plan"));
    ASSERT_TRUE(button->isEnabled());

    button->click();
    ++expected_goal_count;
    ASSERT_TRUE(waitForGui([&]() { return goal_set_server.goalCount() == expected_goal_count; },
                           std::chrono::seconds(2)));
    const auto sent_goal = goal_set_server.lastGoal();
    ASSERT_TRUE(sent_goal.has_value());
    ASSERT_EQ(2u, sent_goal->max_path_lengths.size());
    EXPECT_DOUBLE_EQ(bound, sent_goal->max_path_lengths[0]);
    EXPECT_DOUBLE_EQ(bound, sent_goal->max_path_lengths[1]);
    EXPECT_EQ(0u, single_server.goalCount());

    ASSERT_TRUE(waitForGui(
      [&]() {
        return goal_set_server.terminalSuccessCount() == expected_goal_count &&
               label->text().contains("SUCCEEDED");
      },
      std::chrono::seconds(2)))
      << label->text().toStdString();
    ASSERT_EQ(2, goal_result_table->rowCount());
    EXPECT_DOUBLE_EQ(bound, tableText(*goal_result_table, 0, 3).toDouble());
    EXPECT_DOUBLE_EQ(bound, tableText(*goal_result_table, 1, 3).toDouble());
  }
  EXPECT_EQ(0u, goal_set_server.workerErrorCount());
}

TEST(RaystarPanelInteraction, RejectsGoalSetResultForDifferentMapId) {
  const auto action_name = nextTestName("plan_paths_unused");
  const auto goal_set_action_name = nextTestName("plan_goal_set");
  const auto map_topic = nextTestName("map");
  TestGoalSetActionServer goal_set_server(goal_set_action_name, ServerMode::mismatched_map);
  ASSERT_TRUE(goal_set_server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(goal_set_server.clientNode(), map_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(goal_set_server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(
    context, action_name, map_topic, std::chrono::seconds(2), goal_set_action_name);
  searchModeCombo(*panel)->setCurrentIndex(2);
  auto* label = statusLabel(*panel);
  auto* button = planButton(*panel);
  auto* path_table = resultsTable(*panel);
  auto* goal_table = goalResultsTable(*panel);
  ASSERT_NE(nullptr, label);
  ASSERT_NE(nullptr, button);
  ASSERT_NE(nullptr, path_table);
  ASSERT_NE(nullptr, goal_table);

  const auto map = map_publisher.publishMap();
  ASSERT_TRUE(waitForMapAndPlanReady(*panel, map_publisher, map, std::chrono::seconds(2)));
  button->click();
  ASSERT_TRUE(
    waitForGui([&]() { return label->text().contains(QStringLiteral("different cached map")); },
               std::chrono::seconds(2)))
    << label->text().toStdString();
  EXPECT_EQ(0, path_table->rowCount());
  EXPECT_EQ(0, goal_table->rowCount());
  ASSERT_TRUE(waitForGui([&]() { return goal_set_server.terminalSuccessCount() == 1u; },
                         std::chrono::seconds(2)));
  EXPECT_EQ(0u, goal_set_server.workerErrorCount());
}

TEST(RaystarPanelInteraction, RejectsGoalSetResultWithMismatchedInnerBoundEcho) {
  const auto action_name = nextTestName("plan_paths_unused");
  const auto goal_set_action_name = nextTestName("plan_goal_set");
  const auto map_topic = nextTestName("map");
  TestGoalSetActionServer goal_set_server(goal_set_action_name, ServerMode::inner_bound_mismatch);
  ASSERT_TRUE(goal_set_server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(goal_set_server.clientNode(), map_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(goal_set_server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(
    context, action_name, map_topic, std::chrono::seconds(2), goal_set_action_name);
  searchModeCombo(*panel)->setCurrentIndex(2);
  auto* label = statusLabel(*panel);
  auto* button = planButton(*panel);
  auto* path_table = resultsTable(*panel);
  auto* goal_table = goalResultsTable(*panel);
  ASSERT_NE(nullptr, label);
  ASSERT_NE(nullptr, button);
  ASSERT_NE(nullptr, path_table);
  ASSERT_NE(nullptr, goal_table);

  const auto map = map_publisher.publishMap();
  ASSERT_TRUE(waitForMapAndPlanReady(*panel, map_publisher, map, std::chrono::seconds(2)));
  button->click();
  ASSERT_TRUE(
    waitForGui([&]() { return label->text().startsWith("Failed:"); }, std::chrono::seconds(2)))
    << label->text().toStdString();
  EXPECT_TRUE(label->text().contains("multi-goal"));
  EXPECT_EQ(0, path_table->rowCount());
  EXPECT_EQ(0, goal_table->rowCount());
  ASSERT_TRUE(waitForGui([&]() { return goal_set_server.terminalSuccessCount() == 1u; },
                         std::chrono::seconds(2)));
  EXPECT_EQ(0u, goal_set_server.workerErrorCount());
}

TEST(RaystarPanelInteraction, PartialGoalSetPathsAreShownAsNotPublished) {
  const auto action_name = nextTestName("plan_paths_unused");
  const auto goal_set_action_name = nextTestName("plan_goal_set");
  const auto map_topic = nextTestName("map");
  TestGoalSetActionServer goal_set_server(goal_set_action_name, ServerMode::partial_output);
  ASSERT_TRUE(goal_set_server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(goal_set_server.clientNode(), map_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(goal_set_server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(
    context, action_name, map_topic, std::chrono::seconds(2), goal_set_action_name);
  searchModeCombo(*panel)->setCurrentIndex(2);
  auto* label = statusLabel(*panel);
  auto* button = planButton(*panel);
  auto* path_table = resultsTable(*panel);
  auto* goal_table = goalResultsTable(*panel);
  ASSERT_NE(nullptr, label);
  ASSERT_NE(nullptr, button);
  ASSERT_NE(nullptr, path_table);
  ASSERT_NE(nullptr, goal_table);

  const auto map = map_publisher.publishMap();
  ASSERT_TRUE(waitForMapAndPlanReady(*panel, map_publisher, map, std::chrono::seconds(2)));
  button->click();
  ASSERT_TRUE(waitForGui(
    [&]() {
      return label->text().contains("SUCCEEDED") && label->text().contains("Partial output");
    },
    std::chrono::seconds(2)))
    << label->text().toStdString();

  ASSERT_EQ(1, goal_table->rowCount());
  EXPECT_EQ(QString("Partial output (max paths)"), tableText(*goal_table, 0, 4));
  EXPECT_EQ(QString("no"), tableText(*goal_table, 0, 5));
  EXPECT_EQ(QString("1 / 2"), tableText(*goal_table, 0, 6));
  ASSERT_EQ(1, path_table->rowCount());
  EXPECT_EQ(QString("not published"), tableText(*path_table, 0, 2));
  EXPECT_FALSE(tableText(*path_table, 0, 2).startsWith("path_"));
  EXPECT_TRUE(label->text().contains("Request satisfied: no"));
  EXPECT_TRUE(label->text().contains("path output complete: no"));
  EXPECT_TRUE(label->text().contains("limits: max paths"));
  EXPECT_FALSE(label->text().startsWith("Failed:"));
  EXPECT_FALSE(label->text().contains("does not preserve goal associations"));
  EXPECT_FALSE(label->text().contains("Warning:"));
  ASSERT_TRUE(waitForGui([&]() { return goal_set_server.terminalSuccessCount() == 1u; },
                         std::chrono::seconds(2)));
  EXPECT_EQ(0u, goal_set_server.workerErrorCount());
}

TEST(RaystarPanelInteraction, CapturesStartAndConsecutiveMultiGoalsFromClickedPoint) {
  const auto action_name = nextTestName("plan_paths");
  const auto map_topic = nextTestName("map");
  const auto clicked_point_topic = nextTestName("clicked_point");
  const auto unused_goal_set_action = nextTestName("plan_goal_set_unused");
  TestActionServer server(action_name, ServerMode::success);
  ASSERT_TRUE(server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(server.clientNode(), map_topic);
  TestClickedPointPublisher point_publisher(server.clientNode(), clicked_point_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(context,
                                    action_name,
                                    map_topic,
                                    std::chrono::seconds(2),
                                    unused_goal_set_action,
                                    clicked_point_topic);
  auto* label = statusLabel(*panel);
  auto* capture_start = panel->findChild<QPushButton*>("capture_start_button");
  auto* capture_goal = panel->findChild<QPushButton*>("capture_goal_button");
  auto* start_x = panel->findChild<QLineEdit*>("start_x_edit");
  auto* start_y = panel->findChild<QLineEdit*>("start_y_edit");
  auto* goals = multiGoalTable(*panel);
  ASSERT_NE(nullptr, label);
  ASSERT_NE(nullptr, capture_start);
  ASSERT_NE(nullptr, capture_goal);
  ASSERT_NE(nullptr, start_x);
  ASSERT_NE(nullptr, start_y);
  ASSERT_NE(nullptr, goals);

  const auto map = map_publisher.publishMap();
  ASSERT_TRUE(waitForMapAndPlanReady(*panel, map_publisher, map, std::chrono::seconds(2)));
  ASSERT_TRUE(waitForGui([&]() { return point_publisher.subscriptionCount() > 0u; },
                         std::chrono::seconds(2)));

  capture_start->click();
  ASSERT_TRUE(capture_start->isChecked());
  point_publisher.publish(2.75, -4.125);
  ASSERT_TRUE(waitForGui([&]() { return label->text().contains("Start point captured"); },
                         std::chrono::seconds(2)))
    << label->text().toStdString();
  EXPECT_DOUBLE_EQ(2.75, start_x->text().toDouble());
  EXPECT_DOUBLE_EQ(-4.125, start_y->text().toDouble());
  EXPECT_FALSE(capture_start->isChecked());

  searchModeCombo(*panel)->setCurrentIndex(2);
  goals->setRowCount(0);
  maxPathLengthSpinbox(*panel)->setValue(12.75);
  capture_goal->click();
  ASSERT_TRUE(capture_goal->isChecked());
  point_publisher.publish(6.5, 7.25);
  ASSERT_TRUE(waitForGui([&]() { return goals->rowCount() == 1; }, std::chrono::seconds(2)));
  point_publisher.publish(-1.5, 9.0);
  ASSERT_TRUE(waitForGui([&]() { return goals->rowCount() == 2; }, std::chrono::seconds(2)));
  EXPECT_DOUBLE_EQ(6.5, tableText(*goals, 0, 0).toDouble());
  EXPECT_DOUBLE_EQ(7.25, tableText(*goals, 0, 1).toDouble());
  EXPECT_DOUBLE_EQ(12.75, tableText(*goals, 0, 2).toDouble());
  EXPECT_DOUBLE_EQ(-1.5, tableText(*goals, 1, 0).toDouble());
  EXPECT_DOUBLE_EQ(9.0, tableText(*goals, 1, 1).toDouble());
  EXPECT_DOUBLE_EQ(12.75, tableText(*goals, 1, 2).toDouble());
  EXPECT_TRUE(capture_goal->isChecked());
  capture_goal->click();
  EXPECT_FALSE(capture_goal->isChecked());
  EXPECT_TRUE(label->text().contains("Goal capture stopped"));
}

TEST(RaystarPanelInteraction, MultiGoalDestructionCancelsAndLeavesNoQueuedUiCallback) {
  const auto action_name = nextTestName("plan_paths_unused");
  const auto goal_set_action_name = nextTestName("plan_goal_set");
  const auto map_topic = nextTestName("map");
  TestGoalSetActionServer goal_set_server(
    goal_set_action_name, ServerMode::delayed, std::chrono::seconds(2));
  ASSERT_TRUE(goal_set_server.waitForClient(std::chrono::seconds(2)));
  TestMapPublisher map_publisher(goal_set_server.clientNode(), map_topic);
  auto abstraction = std::make_shared<TestRosNodeAbstraction>(goal_set_server.clientNode());
  TestDisplayContext context(abstraction);
  auto panel = makeInitializedPanel(
    context, action_name, map_topic, std::chrono::seconds(5), goal_set_action_name);
  searchModeCombo(*panel)->setCurrentIndex(2);
  auto* label = statusLabel(*panel);
  auto* button = planButton(*panel);
  ASSERT_NE(nullptr, label);
  ASSERT_NE(nullptr, button);
  const auto map = map_publisher.publishMap();
  ASSERT_TRUE(waitForMapAndPlanReady(*panel, map_publisher, map, std::chrono::seconds(2)));

  button->click();
  ASSERT_TRUE(
    waitForGui([&]() { return goal_set_server.goalCount() > 0u; }, std::chrono::seconds(2)));
  panel.reset();

  ASSERT_TRUE(waitForGui([&]() { return goal_set_server.cancelRequestCount() > 0u; },
                         std::chrono::seconds(2)));
  ASSERT_TRUE(waitForGui([&]() { return goal_set_server.terminalCanceledCount() > 0u; },
                         std::chrono::seconds(2)));
  EXPECT_EQ(0u, goal_set_server.workerErrorCount());
}

}  // namespace

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  rclcpp::init(argc, argv);
  QApplication application(argc, argv);
  const int result = RUN_ALL_TESTS();
  rclcpp::shutdown();
  return result;
}
