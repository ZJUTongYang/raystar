#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <thread>
#include <string>
#include <utility>

#include <QApplication>
#include <QComboBox>
#include <QCoreApplication>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QTableWidget>
#include <QString>
#include <gtest/gtest.h>
#include <geometry_msgs/msg/pose_stamped.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <raystar_interfaces/action/plan_raystar_paths.hpp>
#include <raystar_interfaces/map_identity.hpp>
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

enum class ServerMode { success, mismatched_map, delayed };

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
      [this](const rclcpp_action::GoalUUID&,
             std::shared_ptr<const PlanningAction::Goal> goal) {
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
        workers_.emplace_back([this, goal_handle]() {
          const auto deadline = std::chrono::steady_clock::now() + execution_delay_;
          while (std::chrono::steady_clock::now() < deadline) {
            if (goal_handle->is_canceling()) {
              finishCanceled(goal_handle);
              return;
            }
            std::unique_lock<std::mutex> lock(worker_mutex_);
            worker_cv_.wait_for(lock, std::chrono::milliseconds(5));
          }
          if (goal_handle->is_canceling()) {
            finishCanceled(goal_handle);
            return;
          }
          auto result = std::make_shared<PlanningAction::Result>();
          result->success = mode_ == ServerMode::success;
          result->result_info.status = raystar_interfaces::msg::PlanningResultInfo::STATUS_COMPLETE;
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
        });
      });
    executor_.add_node(node_);
    executor_.add_node(client_node_);
    spin_thread_ = std::thread([this]() { executor_.spin(); });
  }

  ~TestActionServer() {
    executor_.cancel();
    if (spin_thread_.joinable()) {
      spin_thread_.join();
    }
    worker_cv_.notify_all();
    {
      std::lock_guard<std::mutex> lock(worker_mutex_);
      stopping_ = true;
    }
    worker_cv_.notify_all();
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
  std::chrono::milliseconds timeout) {
  auto panel = std::make_unique<raystar_rviz_plugins::RaystarPanel>(nullptr, timeout);
  panel->findChild<QLineEdit*>("action_name_edit")->setText(QString::fromStdString(action_name));
  panel->findChild<QLineEdit*>("map_topic_edit")->setText(QString::fromStdString(map_topic));
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

QComboBox* searchModeCombo(raystar_rviz_plugins::RaystarPanel& panel) {
  return panel.findChild<QComboBox*>("search_mode_combo");
}

QDoubleSpinBox* maxPathLengthSpinbox(raystar_rviz_plugins::RaystarPanel& panel) {
  return panel.findChild<QDoubleSpinBox*>("max_path_length_spinbox");
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
  QString max_path_length;
  int search_mode = -1;
  bool boolean = true;
  ASSERT_TRUE(output.mapGetString("action_name", &action_name));
  EXPECT_EQ(QString("/raystar/plan_paths"), action_name);
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
  ASSERT_TRUE(waitForGui([&]() { return server.lastGoal().has_value(); },
                         std::chrono::seconds(2)));
  const auto goal = server.lastGoal();
  ASSERT_TRUE(goal.has_value());
  EXPECT_EQ(goal->search_mode, PlanningAction::Goal::SEARCH_MODE_ALL_WITHIN_LENGTH);
  EXPECT_EQ(goal->k, 0);
  EXPECT_DOUBLE_EQ(goal->max_path_length, 37.25);
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
