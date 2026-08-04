#pragma once

#include <rviz_common/panel.hpp>
#include <rclcpp/rclcpp.hpp>
#include <rclcpp_action/rclcpp_action.hpp>
#include <geometry_msgs/msg/point_stamped.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <raystar_interfaces/action/plan_raystar_goal_set.hpp>
#include <raystar_interfaces/action/plan_raystar_paths.hpp>
#include <raystar_interfaces/environment_identity.hpp>
#include <raystar_interfaces/map_identity.hpp>

#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QCheckBox>
#include <QLabel>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QTimer>

#include <chrono>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <optional>
#include <string>

namespace raystar_rviz_plugins {

class RaystarPanel : public rviz_common::Panel {
  Q_OBJECT

public:
  explicit RaystarPanel(QWidget* parent = nullptr);
  RaystarPanel(QWidget* parent, std::chrono::milliseconds request_timeout);
  ~RaystarPanel() override;

  void onInitialize() override;
  void save(rviz_common::Config config) const override;
  void load(const rviz_common::Config& config) override;

private Q_SLOTS:
  void onPlanClicked();
  void processCallbacks();
  void onMapTopicChanged(const QString& topic);
  void onClickedPointTopicChanged(const QString& topic);
  void onActionNameChanged(const QString& action_name);
  void onGoalSetActionNameChanged(const QString& action_name);
  void onSearchModeChanged(int index);
  void onMaxPathLengthChanged(double max_path_length);
  void onAddGoalClicked();
  void onRemoveGoalsClicked();
  void onCancelClicked();
  void onCaptureStartClicked(bool checked);
  void onCaptureGoalClicked(bool checked);

private:
  using PlanningAction = raystar_interfaces::action::PlanRaystarPaths;
  using PlanningGoalHandle = rclcpp_action::ClientGoalHandle<PlanningAction>;
  using PlanningActionClient = rclcpp_action::Client<PlanningAction>;
  using GoalSetAction = raystar_interfaces::action::PlanRaystarGoalSet;
  using GoalSetGoalHandle = rclcpp_action::ClientGoalHandle<GoalSetAction>;
  using GoalSetActionClient = rclcpp_action::Client<GoalSetAction>;

  enum class CaptureMode : std::uint8_t { none, start_once, goal_once, append_goals };

  struct CallbackState {
    using Map = nav_msgs::msg::OccupancyGrid;

    struct PendingResponse {
      std::uint64_t request_id{0};
      std::uint64_t map_generation{0};
      rclcpp_action::ResultCode result_code{rclcpp_action::ResultCode::UNKNOWN};
      std::shared_ptr<const PlanningAction::Result> planning_response;
      std::shared_ptr<const GoalSetAction::Result> goal_set_response;
      std::string error;
    };

    struct PendingPoint {
      geometry_msgs::msg::PointStamped point;
      CaptureMode mode{CaptureMode::none};
      std::uint64_t map_generation{0};
    };

    // ROS callbacks only touch this state. Qt widgets and latest_map_ are
    // exclusively accessed by the GUI thread in processCallbacks().
    std::mutex mutex;
    bool alive{true};
    // A subscription id distinguishes callbacks from an old topic after the
    // subscription object has been reset. map_generation changes for every
    // accepted map snapshot, including a topic switch.
    std::uint64_t subscription_id{0};
    std::uint64_t point_subscription_id{0};
    std::uint64_t map_generation{0};
    std::shared_ptr<const Map> pending_map;
    std::uint64_t pending_map_generation{0};
    bool pending_map_update{false};
    std::optional<PendingResponse> pending_response;
    std::deque<PendingPoint> pending_points;
    CaptureMode capture_mode{CaptureMode::none};
    std::uint64_t active_request_id{0};
    std::uint64_t active_map_generation{0};
    PlanningGoalHandle::SharedPtr active_planning_goal;
    GoalSetGoalHandle::SharedPtr active_goal_set_goal;
  };

  void setupUi();
  void subscribeToMap();
  void subscribeToClickedPoint();
  void configureActionClient();
  void configureGoalSetActionClient();
  void updatePlanButtonState();
  void updateGoalRowLabels();
  void addGoalRow(double x, double y, double max_path_length);
  void configurePathResultsTable(bool multi_goal);
  void setCaptureMode(CaptureMode mode);
  std::uint64_t invalidateMapState();
  void cancelActiveRequest();
  void clearResults();
  void planGoalSet();
  void displayResponse(const CallbackState::PendingResponse& result);
  void displayPlanningResponse(const CallbackState::PendingResponse& result);
  void displayGoalSetResponse(const CallbackState::PendingResponse& result);

  QLineEdit* start_x_edit_;
  QLineEdit* start_y_edit_;
  QLineEdit* goal_x_edit_;
  QLineEdit* goal_y_edit_;
  QSpinBox* k_spinbox_;
  QComboBox* search_mode_combo_;
  QLabel* max_path_length_label_;
  QDoubleSpinBox* max_path_length_spinbox_;
  QCheckBox* allow_self_crossing_cb_;
  QCheckBox* allow_unknown_cb_;
  QCheckBox* request_debug_cb_;
  QLineEdit* action_name_edit_;
  QLineEdit* goal_set_action_name_edit_;
  QLineEdit* map_topic_edit_;
  QLineEdit* clicked_point_topic_edit_;
  QGroupBox* single_goal_group_;
  QGroupBox* multi_goal_group_;
  QTableWidget* multi_goal_table_;
  QPushButton* add_goal_button_;
  QPushButton* remove_goals_button_;
  QPushButton* capture_start_button_;
  QPushButton* capture_goal_button_;
  QPushButton* plan_button_;
  QPushButton* cancel_button_;
  QLabel* status_label_;
  QGroupBox* goal_results_group_;
  QTableWidget* goal_results_table_;
  QTableWidget* results_table_;
  QTableWidget* node_table_;
  QLabel* timing_label_;

  rclcpp::Subscription<nav_msgs::msg::OccupancyGrid>::SharedPtr map_sub_;
  rclcpp::Subscription<geometry_msgs::msg::PointStamped>::SharedPtr clicked_point_sub_;
  PlanningActionClient::SharedPtr action_client_;
  GoalSetActionClient::SharedPtr goal_set_action_client_;

  std::shared_ptr<CallbackState> callback_state_;
  QTimer* callback_timer_;
  std::shared_ptr<const nav_msgs::msg::OccupancyGrid> latest_map_;
  std::uint64_t map_generation_{0};
  std::uint64_t next_request_id_{0};
  std::uint64_t active_request_id_{0};
  std::uint64_t active_request_map_generation_{0};
  bool request_in_flight_{false};
  std::optional<std::chrono::steady_clock::time_point> request_deadline_;
  std::string action_client_name_;
  std::string goal_set_action_client_name_;
  std::string subscribed_topic_;
  std::string subscribed_clicked_point_topic_;
  bool ros_initialized_{false};
  std::chrono::milliseconds request_timeout_;
};

}  // namespace raystar_rviz_plugins
