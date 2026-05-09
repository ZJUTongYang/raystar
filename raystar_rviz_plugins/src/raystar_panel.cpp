#include "raystar_rviz_plugins/raystar_panel.h"
#include <pluginlib/class_list_macros.hpp>
#include <rviz_common/display_context.hpp>
#include <QHeaderView>

PLUGINLIB_EXPORT_CLASS(raystar_rviz_plugins::RaystarPanel, rviz_common::Panel)

namespace raystar_rviz_plugins
{

RaystarPanel::RaystarPanel(QWidget* parent)
  : rviz_common::Panel(parent),
    start_x_edit_(new QLineEdit("0.0")),
    start_y_edit_(new QLineEdit("0.0")),
    goal_x_edit_(new QLineEdit("1.0")),
    goal_y_edit_(new QLineEdit("1.0")),
    k_spinbox_(new QSpinBox()),
    allow_self_crossing_cb_(new QCheckBox("Allow Self Crossing")),
    allow_unknown_cb_(new QCheckBox("Allow Unknown")),
    map_topic_edit_(new QLineEdit("/map")),
    plan_button_(new QPushButton("Plan")),
    status_label_(new QLabel("Ready")),
    results_table_(new QTableWidget()),
    timing_label_(new QLabel(""))
{
  k_spinbox_->setRange(1, 100);
  k_spinbox_->setValue(5);
  results_table_->setColumnCount(3);
  results_table_->setHorizontalHeaderLabels({"Path", "Cost", "Waypoints"});
  results_table_->horizontalHeader()->setStretchLastSection(true);
  status_label_->setWordWrap(true);
  timing_label_->setWordWrap(true);
  setupUi();
}

void RaystarPanel::setupUi()
{
  auto* main_layout = new QVBoxLayout;

  auto* map_group = new QGroupBox("Map");
  auto* map_layout = new QHBoxLayout;
  map_layout->addWidget(new QLabel("Topic:"));
  map_layout->addWidget(map_topic_edit_);
  map_group->setLayout(map_layout);
  main_layout->addWidget(map_group);

  auto* start_group = new QGroupBox("Start (world coords)");
  auto* start_layout = new QHBoxLayout;
  start_layout->addWidget(new QLabel("X:"));
  start_layout->addWidget(start_x_edit_);
  start_layout->addWidget(new QLabel("Y:"));
  start_layout->addWidget(start_y_edit_);
  start_group->setLayout(start_layout);
  main_layout->addWidget(start_group);

  auto* goal_group = new QGroupBox("Goal (world coords)");
  auto* goal_layout = new QHBoxLayout;
  goal_layout->addWidget(new QLabel("X:"));
  goal_layout->addWidget(goal_x_edit_);
  goal_layout->addWidget(new QLabel("Y:"));
  goal_layout->addWidget(goal_y_edit_);
  goal_group->setLayout(goal_layout);
  main_layout->addWidget(goal_group);

  auto* param_layout = new QHBoxLayout;
  param_layout->addWidget(new QLabel("K:"));
  param_layout->addWidget(k_spinbox_);
  param_layout->addStretch();
  param_layout->addWidget(allow_self_crossing_cb_);
  param_layout->addWidget(allow_unknown_cb_);
  main_layout->addLayout(param_layout);

  main_layout->addWidget(plan_button_);
  main_layout->addWidget(status_label_);
  main_layout->addWidget(results_table_);
  main_layout->addWidget(timing_label_);
  main_layout->addStretch();

  setLayout(main_layout);

  connect(plan_button_, &QPushButton::clicked, this, &RaystarPanel::onPlanClicked);
  connect(map_topic_edit_, &QLineEdit::editingFinished, this, &RaystarPanel::subscribeToMap);
}

void RaystarPanel::onInitialize()
{
  auto node = getDisplayContext()->getRosNodeAbstraction().lock()->get_raw_node();
  service_client_ = node->create_client<raystar_interfaces::srv::GetRaystarPaths>(
    "/raystar/get_raystar_paths");
  subscribeToMap();
}

void RaystarPanel::subscribeToMap()
{
  auto node = getDisplayContext()->getRosNodeAbstraction().lock()->get_raw_node();
  std::string topic = map_topic_edit_->text().toStdString();
  map_sub_ = node->create_subscription<nav_msgs::msg::OccupancyGrid>(
    topic, rclcpp::QoS(1).transient_local().reliable(),
    [this](nav_msgs::msg::OccupancyGrid::SharedPtr msg) {
      latest_map_ = *msg;
      status_label_->setText(
        QString("Map received: %1x%2 @ %.3f m/cell")
          .arg(msg->info.width).arg(msg->info.height).arg(msg->info.resolution));
    });
}

void RaystarPanel::onPlanClicked()
{
  if (!service_client_->service_is_ready()) {
    status_label_->setText("Error: Service /raystar/get_raystar_paths is not ready");
    return;
  }
  if (latest_map_.data.empty()) {
    status_label_->setText("Error: No map received yet. Check map topic.");
    return;
  }

  auto request = std::make_shared<raystar_interfaces::srv::GetRaystarPaths::Request>();
  request->map = latest_map_;
  request->start.header.frame_id = "map";
  request->start.pose.position.x = start_x_edit_->text().toDouble();
  request->start.pose.position.y = start_y_edit_->text().toDouble();
  request->start.pose.orientation.w = 1.0;
  request->goal.header.frame_id = "map";
  request->goal.pose.position.x = goal_x_edit_->text().toDouble();
  request->goal.pose.position.y = goal_y_edit_->text().toDouble();
  request->goal.pose.orientation.w = 1.0;
  request->k = k_spinbox_->value();
  request->allow_self_crossing = allow_self_crossing_cb_->isChecked();
  request->allow_unknown = allow_unknown_cb_->isChecked();

  status_label_->setText("Planning...");
  plan_button_->setEnabled(false);

  auto future = service_client_->async_send_request(request,
    [this](rclcpp::Client<raystar_interfaces::srv::GetRaystarPaths>::SharedFuture f) {
      auto response = f.get();
      plan_button_->setEnabled(true);

      if (!response->success) {
        status_label_->setText("Failed: " + QString::fromStdString(response->message));
        return;
      }

      status_label_->setText(
        QString("Success: %1 paths found").arg(response->paths.size()));

      results_table_->setRowCount(static_cast<int>(response->paths.size()));
      for (int i = 0; i < static_cast<int>(response->paths.size()); ++i) {
        results_table_->setItem(i, 0, new QTableWidgetItem(QString::number(i + 1)));
        results_table_->setItem(i, 1,
          new QTableWidgetItem(QString::number(response->path_costs[i], 'f', 2)));
        results_table_->setItem(i, 2,
          new QTableWidgetItem(QString::number(response->paths[i].poses.size())));
      }
      results_table_->resizeColumnsToContents();
    });
}

void RaystarPanel::save(rviz_common::Config config) const
{
  rviz_common::Panel::save(config);
  config.mapSetValue("map_topic", map_topic_edit_->text());
  config.mapSetValue("start_x", start_x_edit_->text());
  config.mapSetValue("start_y", start_y_edit_->text());
  config.mapSetValue("goal_x", goal_x_edit_->text());
  config.mapSetValue("goal_y", goal_y_edit_->text());
  config.mapSetValue("k", k_spinbox_->value());
}

void RaystarPanel::load(const rviz_common::Config& config)
{
  rviz_common::Panel::load(config);
  QString val;
  if (config.mapGetString("map_topic", &val)) map_topic_edit_->setText(val);
  if (config.mapGetString("start_x", &val)) start_x_edit_->setText(val);
  if (config.mapGetString("start_y", &val)) start_y_edit_->setText(val);
  if (config.mapGetString("goal_x", &val)) goal_x_edit_->setText(val);
  if (config.mapGetString("goal_y", &val)) goal_y_edit_->setText(val);
  int k_val;
  if (config.mapGetInt("k", &k_val)) k_spinbox_->setValue(k_val);
}

}  // namespace raystar_rviz_plugins
