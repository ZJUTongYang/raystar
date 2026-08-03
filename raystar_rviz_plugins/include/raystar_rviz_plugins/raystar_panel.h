#pragma once

#include <rviz_common/panel.hpp>
#include <rclcpp/rclcpp.hpp>
#include <rclcpp_action/rclcpp_action.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <raystar_interfaces/action/plan_raystar_paths.hpp>
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
  void onActionNameChanged(const QString& action_name);
  void onSearchModeChanged(int index);

private:
  using PlanningAction = raystar_interfaces::action::PlanRaystarPaths;
  using PlanningGoalHandle = rclcpp_action::ClientGoalHandle<PlanningAction>;
  using PlanningActionClient = rclcpp_action::Client<PlanningAction>;

  struct CallbackState {
    using Map = nav_msgs::msg::OccupancyGrid;
    using Response = PlanningAction::Result;

    struct PendingResponse {
      std::uint64_t request_id{0};
      std::uint64_t map_generation{0};
      rclcpp_action::ResultCode result_code{rclcpp_action::ResultCode::UNKNOWN};
      std::shared_ptr<const Response> response;
      std::string error;
    };

    // ROS callbacks only touch this state. Qt widgets and latest_map_ are
    // exclusively accessed by the GUI thread in processCallbacks().
    std::mutex mutex;
    bool alive{true};
    // A subscription id distinguishes callbacks from an old topic after the
    // subscription object has been reset. map_generation changes for every
    // accepted map snapshot, including a topic switch.
    std::uint64_t subscription_id{0};
    std::uint64_t map_generation{0};
    std::shared_ptr<const Map> pending_map;
    std::uint64_t pending_map_generation{0};
    bool pending_map_update{false};
    std::optional<PendingResponse> pending_response;
    std::uint64_t active_request_id{0};
    std::uint64_t active_map_generation{0};
    PlanningGoalHandle::SharedPtr active_goal;
  };

  void setupUi();
  void subscribeToMap();
  void configureActionClient();
  void updatePlanButtonState();
  std::uint64_t invalidateMapState();
  void cancelActiveRequest();
  void clearResults();
  void displayResponse(const CallbackState::PendingResponse& result);

  QLineEdit* start_x_edit_;
  QLineEdit* start_y_edit_;
  QLineEdit* goal_x_edit_;
  QLineEdit* goal_y_edit_;
  QSpinBox* k_spinbox_;
  QComboBox* search_mode_combo_;
  QDoubleSpinBox* max_path_length_spinbox_;
  QCheckBox* allow_self_crossing_cb_;
  QCheckBox* allow_unknown_cb_;
  QCheckBox* request_debug_cb_;
  QLineEdit* action_name_edit_;
  QLineEdit* map_topic_edit_;
  QPushButton* plan_button_;
  QLabel* status_label_;
  QTableWidget* results_table_;
  QTableWidget* node_table_;
  QLabel* timing_label_;

  rclcpp::Subscription<nav_msgs::msg::OccupancyGrid>::SharedPtr map_sub_;
  PlanningActionClient::SharedPtr action_client_;

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
  std::string subscribed_topic_;
  bool ros_initialized_{false};
  std::chrono::milliseconds request_timeout_;
};

}  // namespace raystar_rviz_plugins
