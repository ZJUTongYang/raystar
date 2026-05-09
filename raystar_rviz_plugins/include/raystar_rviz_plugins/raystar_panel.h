#pragma once

#include <rviz_common/panel.hpp>
#include <rclcpp/rclcpp.hpp>
#include <nav_msgs/msg/occupancy_grid.hpp>
#include <raystar_interfaces/srv/get_raystar_paths.hpp>

#include <QLineEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QCheckBox>
#include <QLabel>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>

namespace raystar_rviz_plugins
{

class RaystarPanel : public rviz_common::Panel
{
  Q_OBJECT

public:
  explicit RaystarPanel(QWidget* parent = nullptr);
  ~RaystarPanel() override = default;

  void onInitialize() override;
  void save(rviz_common::Config config) const override;
  void load(const rviz_common::Config& config) override;

private Q_SLOTS:
  void onPlanClicked();

private:
  void setupUi();
  void subscribeToMap();

  QLineEdit* start_x_edit_;
  QLineEdit* start_y_edit_;
  QLineEdit* goal_x_edit_;
  QLineEdit* goal_y_edit_;
  QSpinBox* k_spinbox_;
  QCheckBox* allow_self_crossing_cb_;
  QCheckBox* allow_unknown_cb_;
  QLineEdit* map_topic_edit_;
  QPushButton* plan_button_;
  QLabel* status_label_;
  QTableWidget* results_table_;
  QLabel* timing_label_;

  rclcpp::Subscription<nav_msgs::msg::OccupancyGrid>::SharedPtr map_sub_;
  rclcpp::Client<raystar_interfaces::srv::GetRaystarPaths>::SharedPtr service_client_;

  nav_msgs::msg::OccupancyGrid latest_map_;
};

}  // namespace raystar_rviz_plugins
