#include "raystar_rviz_plugins/raystar_panel.h"
#include <pluginlib/class_list_macros.hpp>
#include <rviz_common/display_context.hpp>
#include <QAbstractItemView>
#include <QGridLayout>
#include <QHeaderView>
#include <QItemSelectionModel>
#include <QSignalBlocker>
#include <QStringList>

#include <algorithm>
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

PLUGINLIB_EXPORT_CLASS(raystar_rviz_plugins::RaystarPanel, rviz_common::Panel)

namespace raystar_rviz_plugins {

namespace {

constexpr std::size_t kMaxUiRows = 10000;
constexpr std::size_t kMaxUiGoals = 1024;
constexpr std::size_t kMaxPendingClickedPoints = 1024;
constexpr std::size_t kMaxUiMessageBytes = 2048;
constexpr auto kDefaultRequestTimeout = std::chrono::seconds(60);
constexpr char kDefaultPlanningActionName[] = "/raystar/plan_paths";
constexpr char kDefaultGoalSetActionName[] = "/raystar/plan_goal_set";
constexpr char kDefaultClickedPointTopic[] = "/clicked_point";
constexpr int kTopKMode = 0;
constexpr int kAllWithinLengthMode = 1;
constexpr int kMultiGoalWithinLengthsMode = 2;

using PlanningAction = raystar_interfaces::action::PlanRaystarPaths;
using PlanningGoalHandle = rclcpp_action::ClientGoalHandle<PlanningAction>;
using PlanningActionClient = rclcpp_action::Client<PlanningAction>;
using GoalSetAction = raystar_interfaces::action::PlanRaystarGoalSet;
using GoalSetGoalHandle = rclcpp_action::ClientGoalHandle<GoalSetAction>;
using GoalSetActionClient = rclcpp_action::Client<GoalSetAction>;
using PlanningResultInfo = raystar_interfaces::msg::PlanningResultInfo;

template <typename ClientPtr, typename GoalHandlePtr>
void cancelGoalBestEffort(const ClientPtr& client, const GoalHandlePtr& goal_handle) noexcept {
  if (!client || !goal_handle) {
    return;
  }
  try {
    (void)client->async_cancel_goal(goal_handle);
  } catch (...) {
    // A terminal goal or a client that is already shutting down no longer
    // needs cancellation. Logical generation guards still reject its result.
  }
}

QString transportStateText(rclcpp_action::ResultCode result_code) {
  switch (result_code) {
  case rclcpp_action::ResultCode::SUCCEEDED:
    return QStringLiteral("SUCCEEDED");
  case rclcpp_action::ResultCode::CANCELED:
    return QStringLiteral("CANCELED");
  case rclcpp_action::ResultCode::ABORTED:
    return QStringLiteral("ABORTED");
  default:
    return QStringLiteral("UNKNOWN");
  }
}

std::string boundedStdString(const std::string& value) {
  if (value.size() <= kMaxUiMessageBytes) {
    return value;
  }
  return value.substr(0, kMaxUiMessageBytes) + " ...[truncated]";
}

QString boundedQString(const std::string& value) {
  return QString::fromStdString(boundedStdString(value));
}

int boundedRowCount(std::size_t count) {
  return static_cast<int>(std::min(count, kMaxUiRows));
}

std::uint64_t nextCounter(std::uint64_t value) {
  return value == std::numeric_limits<std::uint64_t>::max() ? 1 : value + 1;
}

QString planningStatusText(std::uint8_t status) {
  switch (status) {
  case PlanningResultInfo::STATUS_COMPLETE:
    return QStringLiteral("Complete");
  case PlanningResultInfo::STATUS_FEWER_PATHS:
    return QStringLiteral("Complete: fewer paths exist");
  case PlanningResultInfo::STATUS_NO_PATH:
    return QStringLiteral("No path");
  case PlanningResultInfo::STATUS_PARTIAL_SEARCH:
    return QStringLiteral("Partial search");
  case PlanningResultInfo::STATUS_PARTIAL_OUTPUT:
    return QStringLiteral("Partial output");
  case PlanningResultInfo::STATUS_INVALID_REQUEST:
    return QStringLiteral("Invalid request");
  case PlanningResultInfo::STATUS_INVALID_CONFIGURATION:
    return QStringLiteral("Invalid server configuration");
  case PlanningResultInfo::STATUS_BUSY:
    return QStringLiteral("Planner busy");
  case PlanningResultInfo::STATUS_CANCELLED:
    return QStringLiteral("Canceled");
  case PlanningResultInfo::STATUS_FAILED:
    return QStringLiteral("Failed");
  case PlanningResultInfo::STATUS_UNKNOWN:
  default:
    return QStringLiteral("Unknown status");
  }
}

QString planningLimitsText(std::uint16_t limits) {
  if (limits == PlanningResultInfo::LIMIT_NONE) {
    return QStringLiteral("none");
  }
  QStringList names;
  if ((limits & PlanningResultInfo::LIMIT_TIMEOUT) != 0u) {
    names.push_back(QStringLiteral("timeout"));
  }
  if ((limits & PlanningResultInfo::LIMIT_MAX_NODES) != 0u) {
    names.push_back(QStringLiteral("max nodes"));
  }
  if ((limits & PlanningResultInfo::LIMIT_MAX_PATH_POINTS) != 0u) {
    names.push_back(QStringLiteral("max path points"));
  }
  if ((limits & PlanningResultInfo::LIMIT_MAX_RESPONSE_BYTES) != 0u) {
    names.push_back(QStringLiteral("max response bytes"));
  }
  if ((limits & PlanningResultInfo::LIMIT_CANCELLED) != 0u) {
    names.push_back(QStringLiteral("canceled"));
  }
  if ((limits & PlanningResultInfo::LIMIT_MAX_PATHS) != 0u) {
    names.push_back(QStringLiteral("max paths"));
  }
  return names.isEmpty() ? QStringLiteral("unknown") : names.join(", ");
}

}  // namespace

RaystarPanel::RaystarPanel(QWidget* parent) : RaystarPanel(parent, kDefaultRequestTimeout) {}

RaystarPanel::RaystarPanel(QWidget* parent, std::chrono::milliseconds request_timeout)
  : rviz_common::Panel(parent)
  , start_x_edit_(new QLineEdit("0.0"))
  , start_y_edit_(new QLineEdit("0.0"))
  , goal_x_edit_(new QLineEdit("1.0"))
  , goal_y_edit_(new QLineEdit("1.0"))
  , k_spinbox_(new QSpinBox())
  , search_mode_combo_(new QComboBox())
  , max_path_length_label_(new QLabel("Max length (m):"))
  , max_path_length_spinbox_(new QDoubleSpinBox())
  , allow_self_crossing_cb_(new QCheckBox("Allow Self Crossing"))
  , allow_unknown_cb_(new QCheckBox("Allow Unknown"))
  , request_debug_cb_(new QCheckBox("Request Debug"))
  , action_name_edit_(new QLineEdit(kDefaultPlanningActionName))
  , goal_set_action_name_edit_(new QLineEdit(kDefaultGoalSetActionName))
  , map_topic_edit_(new QLineEdit("/map"))
  , clicked_point_topic_edit_(new QLineEdit(kDefaultClickedPointTopic))
  , single_goal_group_(new QGroupBox("Single goal (world coordinates)"))
  , multi_goal_group_(new QGroupBox("Multi-goal targets and inclusive budgets"))
  , multi_goal_table_(new QTableWidget())
  , add_goal_button_(new QPushButton("Add goal"))
  , remove_goals_button_(new QPushButton("Remove selected"))
  , capture_start_button_(new QPushButton("Capture start"))
  , capture_goal_button_(new QPushButton("Capture goal"))
  , plan_button_(new QPushButton("Plan"))
  , cancel_button_(new QPushButton("Cancel"))
  , status_label_(new QLabel("Ready"))
  , goal_results_group_(new QGroupBox("Per-goal completion"))
  , goal_results_table_(new QTableWidget())
  , results_table_(new QTableWidget())
  , node_table_(new QTableWidget())
  , timing_label_(new QLabel(""))
  , callback_state_(std::make_shared<CallbackState>())
  , callback_timer_(new QTimer(this))
  , request_timeout_(request_timeout > std::chrono::milliseconds::zero() ? request_timeout
                                                                         : kDefaultRequestTimeout) {
  k_spinbox_->setRange(1, 100);
  k_spinbox_->setValue(5);
  search_mode_combo_->addItem("Single: Top K");
  search_mode_combo_->addItem("Single: All within length");
  search_mode_combo_->addItem("Multi-goal: All within lengths");
  max_path_length_spinbox_->setRange(1.0e-9, 1.0e9);
  max_path_length_spinbox_->setDecimals(9);
  max_path_length_spinbox_->setValue(10.0);
  max_path_length_spinbox_->setEnabled(false);
  max_path_length_spinbox_->setToolTip(
    "Single-goal length budget. In multi-goal mode, changing this value sets "
    "every existing goal budget; edit table cells afterward for independent budgets.");
  start_x_edit_->setObjectName("start_x_edit");
  start_y_edit_->setObjectName("start_y_edit");
  goal_x_edit_->setObjectName("goal_x_edit");
  goal_y_edit_->setObjectName("goal_y_edit");
  k_spinbox_->setObjectName("k_spinbox");
  search_mode_combo_->setObjectName("search_mode_combo");
  max_path_length_label_->setObjectName("max_path_length_label");
  max_path_length_spinbox_->setObjectName("max_path_length_spinbox");
  allow_self_crossing_cb_->setObjectName("allow_self_crossing_checkbox");
  allow_unknown_cb_->setObjectName("allow_unknown_checkbox");
  request_debug_cb_->setObjectName("request_debug_checkbox");
  action_name_edit_->setObjectName("action_name_edit");
  goal_set_action_name_edit_->setObjectName("goal_set_action_name_edit");
  map_topic_edit_->setObjectName("map_topic_edit");
  clicked_point_topic_edit_->setObjectName("clicked_point_topic_edit");
  multi_goal_group_->setObjectName("multi_goal_group");
  multi_goal_table_->setObjectName("multi_goal_table");
  add_goal_button_->setObjectName("add_goal_button");
  remove_goals_button_->setObjectName("remove_goals_button");
  capture_start_button_->setObjectName("capture_start_button");
  capture_goal_button_->setObjectName("capture_goal_button");
  cancel_button_->setObjectName("cancel_button");
  goal_results_group_->setObjectName("goal_results_group");
  goal_results_table_->setObjectName("goal_results_table");
  capture_start_button_->setCheckable(true);
  capture_goal_button_->setCheckable(true);
  multi_goal_table_->setColumnCount(3);
  multi_goal_table_->setHorizontalHeaderLabels({"X", "Y", "Budget (m)"});
  multi_goal_table_->setToolTip(
    "Inclusive per-goal budgets. Use Set all budgets for a bulk update, then edit "
    "individual cells when goals need different limits.");
  multi_goal_table_->horizontalHeader()->setStretchLastSection(true);
  multi_goal_table_->setSelectionBehavior(QAbstractItemView::SelectRows);
  multi_goal_table_->setSelectionMode(QAbstractItemView::ExtendedSelection);
  goal_results_table_->setColumnCount(7);
  goal_results_table_->setHorizontalHeaderLabels(
    {"Goal", "X", "Y", "Budget (m)", "Status", "Complete", "Paths"});
  goal_results_table_->horizontalHeader()->setStretchLastSection(true);
  configurePathResultsTable(false);
  results_table_->horizontalHeader()->setStretchLastSection(true);
  node_table_->setColumnCount(3);
  node_table_->setHorizontalHeaderLabels({"Node", "G-cost (m)", "F-cost (m)"});
  node_table_->horizontalHeader()->setStretchLastSection(true);
  plan_button_->setEnabled(false);
  cancel_button_->setEnabled(false);
  plan_button_->setObjectName("plan_button");
  status_label_->setObjectName("status_label");
  results_table_->setObjectName("results_table");
  node_table_->setObjectName("node_table");
  timing_label_->setObjectName("timing_label");
  status_label_->setWordWrap(true);
  timing_label_->setWordWrap(true);
  setupUi();
  addGoalRow(1.0, 1.0, 10.0);
  onSearchModeChanged(search_mode_combo_->currentIndex());

  // ROS callbacks run on RViz's executor thread.  Keep all Qt access on the
  // GUI thread by draining a small, mutex-protected handoff periodically.
  callback_timer_->setInterval(50);
  connect(callback_timer_, &QTimer::timeout, this, &RaystarPanel::processCallbacks);
  callback_timer_->start();
}

RaystarPanel::~RaystarPanel() {
  // Stop GUI polling before invalidating the handoff.  ROS callbacks do not
  // capture `this`; they can safely finish while the shared state is alive.
  if (callback_timer_) {
    callback_timer_->stop();
  }
  if (callback_state_) {
    std::lock_guard<std::mutex> lock(callback_state_->mutex);
    callback_state_->alive = false;
    callback_state_->pending_map.reset();
    callback_state_->pending_map_update = false;
    callback_state_->pending_response.reset();
    callback_state_->pending_points.clear();
    callback_state_->capture_mode = CaptureMode::none;
    callback_state_->active_request_id = 0;
  }
  cancelActiveRequest();
  map_sub_.reset();
  clicked_point_sub_.reset();
  action_client_.reset();
  goal_set_action_client_.reset();
  callback_state_.reset();
}

void RaystarPanel::setupUi() {
  auto* main_layout = new QVBoxLayout;

  auto* action_group = new QGroupBox("Planner Actions");
  auto* action_layout = new QGridLayout;
  action_layout->addWidget(new QLabel("Single goal:"), 0, 0);
  action_layout->addWidget(action_name_edit_, 0, 1);
  action_layout->addWidget(new QLabel("Multi-goal:"), 1, 0);
  action_layout->addWidget(goal_set_action_name_edit_, 1, 1);
  action_group->setLayout(action_layout);
  main_layout->addWidget(action_group);

  auto* map_group = new QGroupBox("Map");
  auto* map_layout = new QGridLayout;
  map_layout->addWidget(new QLabel("OccupancyGrid:"), 0, 0);
  map_layout->addWidget(map_topic_edit_, 0, 1);
  map_layout->addWidget(new QLabel("Clicked point:"), 1, 0);
  map_layout->addWidget(clicked_point_topic_edit_, 1, 1);
  map_group->setLayout(map_layout);
  main_layout->addWidget(map_group);

  auto* start_group = new QGroupBox("Start (world coords)");
  auto* start_layout = new QHBoxLayout;
  start_layout->addWidget(new QLabel("X:"));
  start_layout->addWidget(start_x_edit_);
  start_layout->addWidget(new QLabel("Y:"));
  start_layout->addWidget(start_y_edit_);
  start_layout->addWidget(capture_start_button_);
  start_group->setLayout(start_layout);
  main_layout->addWidget(start_group);

  auto* goal_layout = new QHBoxLayout;
  goal_layout->addWidget(new QLabel("X:"));
  goal_layout->addWidget(goal_x_edit_);
  goal_layout->addWidget(new QLabel("Y:"));
  goal_layout->addWidget(goal_y_edit_);
  single_goal_group_->setLayout(goal_layout);
  main_layout->addWidget(single_goal_group_);

  auto* multi_goal_layout = new QVBoxLayout;
  multi_goal_layout->addWidget(multi_goal_table_);
  auto* multi_goal_buttons = new QHBoxLayout;
  multi_goal_buttons->addWidget(add_goal_button_);
  multi_goal_buttons->addWidget(remove_goals_button_);
  multi_goal_buttons->addStretch();
  multi_goal_layout->addLayout(multi_goal_buttons);
  multi_goal_group_->setLayout(multi_goal_layout);
  main_layout->addWidget(multi_goal_group_);

  auto* goal_capture_layout = new QHBoxLayout;
  goal_capture_layout->addWidget(new QLabel("Publish Point tool:"));
  goal_capture_layout->addWidget(capture_goal_button_);
  goal_capture_layout->addStretch();
  main_layout->addLayout(goal_capture_layout);

  auto* param_layout = new QHBoxLayout;
  param_layout->addWidget(new QLabel("Mode:"));
  param_layout->addWidget(search_mode_combo_);
  param_layout->addWidget(new QLabel("K:"));
  param_layout->addWidget(k_spinbox_);
  param_layout->addWidget(max_path_length_label_);
  param_layout->addWidget(max_path_length_spinbox_);
  param_layout->addStretch();
  param_layout->addWidget(allow_self_crossing_cb_);
  param_layout->addWidget(allow_unknown_cb_);
  param_layout->addWidget(request_debug_cb_);
  main_layout->addLayout(param_layout);

  auto* request_buttons = new QHBoxLayout;
  request_buttons->addWidget(plan_button_);
  request_buttons->addWidget(cancel_button_);
  main_layout->addLayout(request_buttons);
  main_layout->addWidget(status_label_);
  auto* goal_results_layout = new QVBoxLayout;
  goal_results_layout->addWidget(goal_results_table_);
  goal_results_group_->setLayout(goal_results_layout);
  main_layout->addWidget(goal_results_group_);
  main_layout->addWidget(results_table_);
  main_layout->addWidget(new QLabel("Tree nodes:"));
  main_layout->addWidget(node_table_);
  main_layout->addWidget(timing_label_);
  main_layout->addStretch();

  setLayout(main_layout);

  connect(plan_button_, &QPushButton::clicked, this, &RaystarPanel::onPlanClicked);
  connect(map_topic_edit_, &QLineEdit::textChanged, this, &RaystarPanel::onMapTopicChanged);
  connect(map_topic_edit_, &QLineEdit::editingFinished, this, &RaystarPanel::subscribeToMap);
  connect(clicked_point_topic_edit_,
          &QLineEdit::textChanged,
          this,
          &RaystarPanel::onClickedPointTopicChanged);
  connect(clicked_point_topic_edit_,
          &QLineEdit::editingFinished,
          this,
          &RaystarPanel::subscribeToClickedPoint);
  connect(action_name_edit_, &QLineEdit::textChanged, this, &RaystarPanel::onActionNameChanged);
  connect(goal_set_action_name_edit_,
          &QLineEdit::textChanged,
          this,
          &RaystarPanel::onGoalSetActionNameChanged);
  connect(search_mode_combo_,
          qOverload<int>(&QComboBox::currentIndexChanged),
          this,
          &RaystarPanel::onSearchModeChanged);
  connect(max_path_length_spinbox_,
          qOverload<double>(&QDoubleSpinBox::valueChanged),
          this,
          &RaystarPanel::onMaxPathLengthChanged);
  connect(
    action_name_edit_, &QLineEdit::editingFinished, this, &RaystarPanel::configureActionClient);
  connect(goal_set_action_name_edit_,
          &QLineEdit::editingFinished,
          this,
          &RaystarPanel::configureGoalSetActionClient);
  connect(add_goal_button_, &QPushButton::clicked, this, &RaystarPanel::onAddGoalClicked);
  connect(remove_goals_button_, &QPushButton::clicked, this, &RaystarPanel::onRemoveGoalsClicked);
  connect(cancel_button_, &QPushButton::clicked, this, &RaystarPanel::onCancelClicked);
  connect(capture_start_button_, &QPushButton::toggled, this, &RaystarPanel::onCaptureStartClicked);
  connect(capture_goal_button_, &QPushButton::toggled, this, &RaystarPanel::onCaptureGoalClicked);
}

void RaystarPanel::onSearchModeChanged(int index) {
  const bool multi_goal = index == kMultiGoalWithinLengthsMode;
  const bool bounded = index != kTopKMode;
  if (request_in_flight_) {
    cancelActiveRequest();
    status_label_->setText("Planning mode changed; the previous goal was canceled.");
  }
  setCaptureMode(CaptureMode::none);
  k_spinbox_->setEnabled(!bounded);
  max_path_length_spinbox_->setEnabled(bounded);
  max_path_length_label_->setText(multi_goal ? "Set all budgets (m):" : "Max length (m):");
  single_goal_group_->setVisible(!multi_goal);
  multi_goal_group_->setVisible(multi_goal);
  goal_results_group_->setVisible(multi_goal);
  capture_goal_button_->setText(multi_goal ? "Capture goals" : "Capture goal");
  clearResults();
  updatePlanButtonState();
}

void RaystarPanel::onMaxPathLengthChanged(double max_path_length) {
  const int mode = search_mode_combo_->currentIndex();
  if (mode == kTopKMode || !std::isfinite(max_path_length) || max_path_length <= 0.0) {
    return;
  }

  const bool canceled_request = request_in_flight_;
  if (canceled_request) {
    cancelActiveRequest();
  }

  if (mode == kMultiGoalWithinLengthsMode) {
    const QSignalBlocker blocker(multi_goal_table_);
    const QString budget_text = QString::number(max_path_length, 'g', 17);
    for (int row = 0; row < multi_goal_table_->rowCount(); ++row) {
      auto* budget_item = multi_goal_table_->item(row, 2);
      if (!budget_item) {
        budget_item = new QTableWidgetItem;
        multi_goal_table_->setItem(row, 2, budget_item);
      }
      budget_item->setText(budget_text);
    }
    clearResults();
    status_label_->setText(
      QString("%1Set all %2 goal budgets to %3 m. Press Plan to refresh paths.")
        .arg(canceled_request ? QStringLiteral("Previous request canceled. ") : QString())
        .arg(multi_goal_table_->rowCount())
        .arg(max_path_length, 0, 'g', 10));
  } else {
    clearResults();
    status_label_->setText(
      QString("%1Maximum path length set to %2 m. Press Plan to refresh paths.")
        .arg(canceled_request ? QStringLiteral("Previous request canceled. ") : QString())
        .arg(max_path_length, 0, 'g', 10));
  }
  updatePlanButtonState();
}

void RaystarPanel::onInitialize() {
  if (ros_initialized_) {
    return;
  }
  auto* context = getDisplayContext();
  if (!context) {
    status_label_->setText("Error: RViz ROS context is unavailable.");
    return;
  }
  auto abstraction = context->getRosNodeAbstraction().lock();
  if (!abstraction) {
    status_label_->setText("Error: RViz ROS node is unavailable.");
    return;
  }
  auto node = abstraction->get_raw_node();
  if (!node) {
    status_label_->setText("Error: RViz ROS node is unavailable.");
    return;
  }
  ros_initialized_ = true;
  subscribeToMap();
  subscribeToClickedPoint();
  configureActionClient();
  configureGoalSetActionClient();
}

void RaystarPanel::updatePlanButtonState() {
  const bool multi_goal = search_mode_combo_->currentIndex() == kMultiGoalWithinLengthsMode;
  const bool action_is_current =
    multi_goal ? goal_set_action_client_ &&
                   goal_set_action_name_edit_->text().toStdString() == goal_set_action_client_name_
               : action_client_ && action_name_edit_->text().toStdString() == action_client_name_;
  const bool map_is_current = latest_map_ && !latest_map_->data.empty() &&
                              map_topic_edit_->text().toStdString() == subscribed_topic_;
  const bool goals_available = !multi_goal || multi_goal_table_->rowCount() > 0;
  plan_button_->setEnabled(ros_initialized_ && action_is_current && map_is_current &&
                           goals_available && !request_in_flight_);
  cancel_button_->setEnabled(request_in_flight_);
}

void RaystarPanel::onActionNameChanged(const QString& action_name) {
  if (action_name.toStdString() == action_client_name_) {
    return;
  }

  // Keep the old client alive until cancellation has been requested.  A goal
  // whose acceptance response arrives later is also canceled by the callback,
  // which captured that same client when the request was sent.
  cancelActiveRequest();
  action_client_.reset();
  action_client_name_.clear();
  clearResults();
  updatePlanButtonState();
  status_label_->setText(ros_initialized_
                           ? "Action name changed; press Enter or leave the field to reconnect."
                           : "Action selected; waiting for RViz ROS initialization.");
}

void RaystarPanel::onGoalSetActionNameChanged(const QString& action_name) {
  if (action_name.toStdString() == goal_set_action_client_name_) {
    return;
  }
  cancelActiveRequest();
  goal_set_action_client_.reset();
  goal_set_action_client_name_.clear();
  clearResults();
  updatePlanButtonState();
  status_label_->setText(
    ros_initialized_
      ? "Multi-goal Action name changed; press Enter or leave the field to reconnect."
      : "Multi-goal Action selected; waiting for RViz ROS initialization.");
}

void RaystarPanel::configureActionClient() {
  if (!ros_initialized_) {
    status_label_->setText("Action selected; waiting for RViz ROS initialization.");
    return;
  }

  const std::string action_name = action_name_edit_->text().toStdString();
  if (action_name.empty()) {
    cancelActiveRequest();
    action_client_.reset();
    action_client_name_.clear();
    updatePlanButtonState();
    status_label_->setText("Error: Action name cannot be empty.");
    return;
  }
  if (action_client_ && action_name == action_client_name_) {
    updatePlanButtonState();
    return;
  }

  auto* context = getDisplayContext();
  auto abstraction = context ? context->getRosNodeAbstraction().lock() : nullptr;
  auto node = abstraction ? abstraction->get_raw_node() : nullptr;
  if (!node) {
    cancelActiveRequest();
    action_client_.reset();
    action_client_name_.clear();
    updatePlanButtonState();
    status_label_->setText("Error: RViz ROS node is unavailable.");
    return;
  }

  cancelActiveRequest();
  action_client_.reset();
  action_client_name_.clear();
  try {
    action_client_ = rclcpp_action::create_client<PlanningAction>(node, action_name);
    action_client_name_ = action_name;
    status_label_->setText(QString("Using planning Action %1").arg(boundedQString(action_name)));
  } catch (const std::exception& exception) {
    action_client_.reset();
    status_label_->setText(
      boundedQString(std::string("Could not create Action client: ") + exception.what()));
  } catch (...) {
    action_client_.reset();
    status_label_->setText("Error: Could not create the Raystar Action client.");
  }
  updatePlanButtonState();
}

void RaystarPanel::configureGoalSetActionClient() {
  if (!ros_initialized_) {
    status_label_->setText("Multi-goal Action selected; waiting for RViz ROS initialization.");
    return;
  }

  const std::string action_name = goal_set_action_name_edit_->text().toStdString();
  if (action_name.empty()) {
    cancelActiveRequest();
    goal_set_action_client_.reset();
    goal_set_action_client_name_.clear();
    updatePlanButtonState();
    status_label_->setText("Error: Multi-goal Action name cannot be empty.");
    return;
  }
  if (goal_set_action_client_ && action_name == goal_set_action_client_name_) {
    updatePlanButtonState();
    return;
  }

  auto* context = getDisplayContext();
  auto abstraction = context ? context->getRosNodeAbstraction().lock() : nullptr;
  auto node = abstraction ? abstraction->get_raw_node() : nullptr;
  if (!node) {
    cancelActiveRequest();
    goal_set_action_client_.reset();
    goal_set_action_client_name_.clear();
    updatePlanButtonState();
    status_label_->setText("Error: RViz ROS node is unavailable.");
    return;
  }

  cancelActiveRequest();
  goal_set_action_client_.reset();
  goal_set_action_client_name_.clear();
  try {
    goal_set_action_client_ = rclcpp_action::create_client<GoalSetAction>(node, action_name);
    goal_set_action_client_name_ = action_name;
    status_label_->setText(QString("Using multi-goal Action %1").arg(boundedQString(action_name)));
  } catch (const std::exception& exception) {
    goal_set_action_client_.reset();
    status_label_->setText(boundedQString(
      std::string("Could not create multi-goal Action client: ") + exception.what()));
  } catch (...) {
    goal_set_action_client_.reset();
    status_label_->setText("Error: Could not create the multi-goal Action client.");
  }
  updatePlanButtonState();
}

std::uint64_t RaystarPanel::invalidateMapState() {
  map_generation_ = nextCounter(map_generation_);
  cancelActiveRequest();
  setCaptureMode(CaptureMode::none);
  latest_map_.reset();
  subscribed_topic_.clear();
  plan_button_->setEnabled(false);
  clearResults();

  std::uint64_t subscription_id = 0;
  if (callback_state_) {
    std::lock_guard<std::mutex> lock(callback_state_->mutex);
    if (!callback_state_->alive) {
      return 0;
    }
    callback_state_->subscription_id = nextCounter(callback_state_->subscription_id);
    subscription_id = callback_state_->subscription_id;
    callback_state_->map_generation = map_generation_;
    callback_state_->pending_map.reset();
    callback_state_->pending_map_update = false;
    callback_state_->pending_response.reset();
    callback_state_->pending_points.clear();
    callback_state_->active_request_id = 0;
    callback_state_->active_map_generation = map_generation_;
  }

  map_sub_.reset();
  return subscription_id;
}

void RaystarPanel::onMapTopicChanged(const QString& topic) {
  if (topic.toStdString() == subscribed_topic_) {
    return;
  }
  invalidateMapState();
  status_label_->setText(ros_initialized_
                           ? "Map topic changed; press Enter or leave the field to subscribe."
                           : "Map topic selected; waiting for RViz ROS initialization.");
}

void RaystarPanel::subscribeToMap() {
  const std::string topic = map_topic_edit_->text().toStdString();
  if (topic.empty()) {
    status_label_->setText("Error: Map topic cannot be empty.");
    return;
  }
  if (map_sub_ && topic == subscribed_topic_) {
    updatePlanButtonState();
    return;
  }

  // Topic changes invalidate both the displayed map and any response for a
  // previous map.  The generation is checked again when the GUI drains the
  // callback handoff, so a late callback cannot resurrect stale state.
  const std::uint64_t subscription_id = invalidateMapState();
  if (subscription_id == 0) {
    return;
  }

  auto* context = getDisplayContext();
  if (!ros_initialized_ || !context) {
    status_label_->setText("Map topic selected; waiting for RViz ROS initialization.");
    return;
  }
  auto abstraction = context->getRosNodeAbstraction().lock();
  if (!abstraction) {
    status_label_->setText("Error: RViz ROS node is unavailable.");
    return;
  }
  auto node = abstraction->get_raw_node();
  if (!node) {
    status_label_->setText("Error: RViz ROS node is unavailable.");
    return;
  }
  const auto weak_state = std::weak_ptr<CallbackState>(callback_state_);
  try {
    map_sub_ = node->create_subscription<nav_msgs::msg::OccupancyGrid>(
      topic,
      rclcpp::QoS(1).transient_local().reliable(),
      [weak_state, subscription_id](nav_msgs::msg::OccupancyGrid::ConstSharedPtr msg) {
        if (!msg) {
          return;
        }
        auto state = weak_state.lock();
        if (!state) {
          return;
        }
        std::shared_ptr<const nav_msgs::msg::OccupancyGrid> snapshot = msg;
        std::lock_guard<std::mutex> lock(state->mutex);
        if (!state->alive || state->subscription_id != subscription_id) {
          return;
        }
        // QoS depth is one, and this handoff deliberately keeps only the
        // newest snapshot to avoid an unbounded queue while RViz is busy.
        state->map_generation = nextCounter(state->map_generation);
        state->pending_map = std::move(snapshot);
        state->pending_map_generation = state->map_generation;
        state->pending_map_update = true;
        // A new map invalidates a request made against the previous map.
        state->active_request_id = 0;
        state->active_map_generation = state->map_generation;
        state->pending_response.reset();
        state->pending_points.clear();
        state->capture_mode = CaptureMode::none;
      });
    subscribed_topic_ = topic;
    status_label_->setText(QString("Waiting for map on %1").arg(boundedQString(topic)));
  } catch (const std::exception&) {
    map_sub_.reset();
    status_label_->setText("Error: Could not subscribe to the selected map topic.");
  }
}

void RaystarPanel::onClickedPointTopicChanged(const QString& topic) {
  if (topic.toStdString() == subscribed_clicked_point_topic_) {
    return;
  }
  setCaptureMode(CaptureMode::none);
  clicked_point_sub_.reset();
  subscribed_clicked_point_topic_.clear();
  if (callback_state_) {
    std::lock_guard<std::mutex> lock(callback_state_->mutex);
    callback_state_->point_subscription_id = nextCounter(callback_state_->point_subscription_id);
    callback_state_->pending_points.clear();
  }
  status_label_->setText(
    ros_initialized_ ? "Clicked-point topic changed; press Enter or leave the field to subscribe."
                     : "Clicked-point topic selected; waiting for RViz ROS initialization.");
}

void RaystarPanel::subscribeToClickedPoint() {
  const std::string topic = clicked_point_topic_edit_->text().toStdString();
  if (topic.empty()) {
    status_label_->setText("Error: Clicked-point topic cannot be empty.");
    return;
  }
  if (clicked_point_sub_ && topic == subscribed_clicked_point_topic_) {
    return;
  }

  auto* context = getDisplayContext();
  auto abstraction = context ? context->getRosNodeAbstraction().lock() : nullptr;
  auto node = abstraction ? abstraction->get_raw_node() : nullptr;
  if (!ros_initialized_ || !node) {
    status_label_->setText("Clicked-point topic selected; waiting for RViz ROS initialization.");
    return;
  }

  setCaptureMode(CaptureMode::none);
  clicked_point_sub_.reset();
  subscribed_clicked_point_topic_.clear();
  std::uint64_t subscription_id = 0;
  if (callback_state_) {
    std::lock_guard<std::mutex> lock(callback_state_->mutex);
    if (!callback_state_->alive) {
      return;
    }
    callback_state_->point_subscription_id = nextCounter(callback_state_->point_subscription_id);
    subscription_id = callback_state_->point_subscription_id;
    callback_state_->pending_points.clear();
  }
  if (subscription_id == 0) {
    return;
  }

  const auto weak_state = std::weak_ptr<CallbackState>(callback_state_);
  try {
    clicked_point_sub_ = node->create_subscription<geometry_msgs::msg::PointStamped>(
      topic,
      rclcpp::QoS(10).reliable(),
      [weak_state, subscription_id](geometry_msgs::msg::PointStamped::ConstSharedPtr msg) {
        if (!msg) {
          return;
        }
        auto state = weak_state.lock();
        if (!state) {
          return;
        }
        std::lock_guard<std::mutex> lock(state->mutex);
        if (!state->alive || state->point_subscription_id != subscription_id ||
            state->capture_mode == CaptureMode::none ||
            state->pending_points.size() >= kMaxPendingClickedPoints) {
          return;
        }
        CallbackState::PendingPoint pending;
        pending.point = *msg;
        pending.mode = state->capture_mode;
        pending.map_generation = state->map_generation;
        state->pending_points.emplace_back(std::move(pending));
        if (state->capture_mode == CaptureMode::start_once ||
            state->capture_mode == CaptureMode::goal_once) {
          state->capture_mode = CaptureMode::none;
        }
      });
    subscribed_clicked_point_topic_ = topic;
  } catch (const std::exception& exception) {
    clicked_point_sub_.reset();
    status_label_->setText(
      boundedQString(std::string("Could not subscribe to clicked points: ") + exception.what()));
  } catch (...) {
    clicked_point_sub_.reset();
    status_label_->setText("Error: Could not subscribe to the clicked-point topic.");
  }
}

void RaystarPanel::updateGoalRowLabels() {
  for (int row = 0; row < multi_goal_table_->rowCount(); ++row) {
    multi_goal_table_->setVerticalHeaderItem(
      row, new QTableWidgetItem(QString::number(static_cast<qlonglong>(row))));
  }
}

void RaystarPanel::addGoalRow(double x, double y, double max_path_length) {
  if (static_cast<std::size_t>(multi_goal_table_->rowCount()) >= kMaxUiGoals) {
    status_label_->setText(QString("Error: The Panel accepts at most %1 goals.")
                             .arg(static_cast<qulonglong>(kMaxUiGoals)));
    return;
  }
  const int row = multi_goal_table_->rowCount();
  multi_goal_table_->insertRow(row);
  multi_goal_table_->setItem(row, 0, new QTableWidgetItem(QString::number(x, 'g', 17)));
  multi_goal_table_->setItem(row, 1, new QTableWidgetItem(QString::number(y, 'g', 17)));
  multi_goal_table_->setItem(
    row, 2, new QTableWidgetItem(QString::number(max_path_length, 'g', 17)));
  updateGoalRowLabels();
  updatePlanButtonState();
}

void RaystarPanel::onAddGoalClicked() {
  bool x_ok = false;
  bool y_ok = false;
  double x = goal_x_edit_->text().toDouble(&x_ok);
  double y = goal_y_edit_->text().toDouble(&y_ok);
  if ((!x_ok || !y_ok || !std::isfinite(x) || !std::isfinite(y)) &&
      multi_goal_table_->rowCount() > 0) {
    const int last = multi_goal_table_->rowCount() - 1;
    x = multi_goal_table_->item(last, 0)->text().toDouble(&x_ok);
    y = multi_goal_table_->item(last, 1)->text().toDouble(&y_ok);
  }
  if (!x_ok || !y_ok || !std::isfinite(x) || !std::isfinite(y)) {
    x = 0.0;
    y = 0.0;
  }
  addGoalRow(x, y, max_path_length_spinbox_->value());
}

void RaystarPanel::onRemoveGoalsClicked() {
  const auto rows = multi_goal_table_->selectionModel()->selectedRows();
  if (rows.empty()) {
    status_label_->setText("Select one or more goal rows to remove.");
    return;
  }
  std::vector<int> indices;
  indices.reserve(static_cast<std::size_t>(rows.size()));
  for (const auto& index : rows) {
    indices.push_back(index.row());
  }
  std::sort(indices.begin(), indices.end(), std::greater<int>());
  indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
  for (const int row : indices) {
    multi_goal_table_->removeRow(row);
  }
  updateGoalRowLabels();
  clearResults();
  updatePlanButtonState();
}

void RaystarPanel::configurePathResultsTable(bool multi_goal) {
  if (multi_goal) {
    results_table_->setColumnCount(5);
    results_table_->setHorizontalHeaderLabels(
      {"Goal", "Path", "Expected marker", "Cost (m)", "Waypoints"});
  } else {
    results_table_->setColumnCount(3);
    results_table_->setHorizontalHeaderLabels({"Path", "Cost (m)", "Waypoints"});
  }
}

void RaystarPanel::setCaptureMode(CaptureMode mode) {
  if (callback_state_) {
    std::lock_guard<std::mutex> lock(callback_state_->mutex);
    if (callback_state_->alive) {
      callback_state_->capture_mode = mode;
      callback_state_->pending_points.clear();
    }
  }
  const QSignalBlocker start_blocker(capture_start_button_);
  const QSignalBlocker goal_blocker(capture_goal_button_);
  capture_start_button_->setChecked(mode == CaptureMode::start_once);
  capture_goal_button_->setChecked(mode == CaptureMode::goal_once ||
                                   mode == CaptureMode::append_goals);
}

void RaystarPanel::onCaptureStartClicked(bool checked) {
  if (!checked) {
    setCaptureMode(CaptureMode::none);
    return;
  }
  if (!latest_map_) {
    setCaptureMode(CaptureMode::none);
    status_label_->setText("Error: Receive a map before capturing a start point.");
    return;
  }
  if (!clicked_point_sub_ ||
      clicked_point_topic_edit_->text().toStdString() != subscribed_clicked_point_topic_) {
    subscribeToClickedPoint();
  }
  if (!clicked_point_sub_) {
    setCaptureMode(CaptureMode::none);
    return;
  }
  setCaptureMode(CaptureMode::start_once);
  status_label_->setText("Select RViz's Publish Point tool and click the start position.");
}

void RaystarPanel::onCaptureGoalClicked(bool checked) {
  if (!checked) {
    setCaptureMode(CaptureMode::none);
    status_label_->setText("Goal capture stopped.");
    return;
  }
  if (!latest_map_) {
    setCaptureMode(CaptureMode::none);
    status_label_->setText("Error: Receive a map before capturing goals.");
    return;
  }
  if (!clicked_point_sub_ ||
      clicked_point_topic_edit_->text().toStdString() != subscribed_clicked_point_topic_) {
    subscribeToClickedPoint();
  }
  if (!clicked_point_sub_) {
    setCaptureMode(CaptureMode::none);
    return;
  }
  const bool multi_goal = search_mode_combo_->currentIndex() == kMultiGoalWithinLengthsMode;
  setCaptureMode(multi_goal ? CaptureMode::append_goals : CaptureMode::goal_once);
  status_label_->setText(
    multi_goal
      ? "Select RViz's Publish Point tool and click each goal; toggle Capture goals to stop."
      : "Select RViz's Publish Point tool and click the goal position.");
}

void RaystarPanel::onCancelClicked() {
  if (!request_in_flight_) {
    return;
  }
  cancelActiveRequest();
  updatePlanButtonState();
  status_label_->setText("Cancellation requested; the displayed result was invalidated.");
}

void RaystarPanel::onPlanClicked() {
  if (search_mode_combo_->currentIndex() == kMultiGoalWithinLengthsMode) {
    planGoalSet();
    return;
  }
  // Drain a map callback that may have arrived just before the click.  This
  // keeps the request snapshot and the map generation in the GUI thread.
  processCallbacks();

  const std::string selected_topic = map_topic_edit_->text().toStdString();
  if (selected_topic.empty() || selected_topic != subscribed_topic_) {
    subscribeToMap();
    return;
  }

  const std::string selected_action = action_name_edit_->text().toStdString();
  if (selected_action.empty() || selected_action != action_client_name_ || !action_client_) {
    configureActionClient();
    return;
  }
  if (request_in_flight_) {
    status_label_->setText("Planning is already in progress.");
    return;
  }
  if (!action_client_->action_server_is_ready()) {
    status_label_->setText(
      QString("Error: Action %1 is not ready").arg(boundedQString(action_client_name_)));
    return;
  }
  if (!latest_map_ || latest_map_->data.empty()) {
    status_label_->setText("Error: No map received yet. Check map topic.");
    return;
  }
  if (latest_map_->header.frame_id.empty()) {
    status_label_->setText("Error: The received map has an empty frame_id.");
    return;
  }

  bool start_x_ok = false;
  bool start_y_ok = false;
  bool goal_x_ok = false;
  bool goal_y_ok = false;
  const double start_x = start_x_edit_->text().toDouble(&start_x_ok);
  const double start_y = start_y_edit_->text().toDouble(&start_y_ok);
  const double goal_x = goal_x_edit_->text().toDouble(&goal_x_ok);
  const double goal_y = goal_y_edit_->text().toDouble(&goal_y_ok);
  if (!start_x_ok || !start_y_ok || !goal_x_ok || !goal_y_ok || !std::isfinite(start_x) ||
      !std::isfinite(start_y) || !std::isfinite(goal_x) || !std::isfinite(goal_y)) {
    status_label_->setText("Error: Start and goal coordinates must be finite numbers.");
    return;
  }

  auto state = callback_state_;
  if (!state) {
    status_label_->setText("Error: Panel callback state is unavailable.");
    return;
  }

  std::uint64_t request_id = 0;
  std::uint64_t request_generation = 0;
  {
    std::lock_guard<std::mutex> lock(state->mutex);
    if (!state->alive) {
      status_label_->setText("Error: Panel is shutting down.");
      return;
    }
    // A ROS callback can publish a newer snapshot between the timer tick and
    // this click.  Do not send a request using the now-stale GUI snapshot.
    if (state->map_generation != map_generation_ || state->pending_map_update) {
      status_label_->setText("Map changed; wait for the new map before planning.");
      return;
    }
    request_generation = state->map_generation;
  }

  PlanningAction::Goal goal;
  try {
    // Identify the immutable GUI snapshot without copying its potentially
    // large data vector into the Action sample. The server accepts the goal
    // only when its cached map has the same deterministic identity.
    goal.map_id = raystar_interfaces::computeMapId(*latest_map_);
    goal.start.header.frame_id = latest_map_->header.frame_id;
    goal.start.pose.position.x = start_x;
    goal.start.pose.position.y = start_y;
    goal.start.pose.orientation.w = 1.0;
    goal.goal.header.frame_id = latest_map_->header.frame_id;
    goal.goal.pose.position.x = goal_x;
    goal.goal.pose.position.y = goal_y;
    goal.goal.pose.orientation.w = 1.0;
    const bool bounded_mode = search_mode_combo_->currentIndex() == 1;
    goal.search_mode = bounded_mode ? PlanningAction::Goal::SEARCH_MODE_ALL_WITHIN_LENGTH
                                    : PlanningAction::Goal::SEARCH_MODE_TOP_K;
    goal.k = bounded_mode ? 0 : k_spinbox_->value();
    goal.max_path_length = bounded_mode ? max_path_length_spinbox_->value() : 0.0;
    goal.allow_self_crossing = allow_self_crossing_cb_->isChecked();
    goal.allow_unknown = allow_unknown_cb_->isChecked();
    goal.include_debug = request_debug_cb_->isChecked();
  } catch (const std::exception& exception) {
    status_label_->setText(
      boundedQString(std::string("Could not prepare planning request: ") + exception.what()));
    return;
  } catch (...) {
    status_label_->setText("Error: Could not prepare planning request.");
    return;
  }

  // Re-check the generation after hashing the map. A callback may have
  // arrived while the potentially large occupancy vector was being scanned.
  {
    std::lock_guard<std::mutex> lock(state->mutex);
    if (!state->alive || state->map_generation != request_generation || state->pending_map_update) {
      status_label_->setText("Map changed; wait for the new map before planning.");
      return;
    }
    next_request_id_ = nextCounter(next_request_id_);
    request_id = next_request_id_;
    state->active_request_id = request_id;
    state->active_map_generation = request_generation;
    state->active_planning_goal.reset();
    state->active_goal_set_goal.reset();
    state->pending_response.reset();
  }

  clearResults();
  status_label_->setText("Planning...");
  plan_button_->setEnabled(false);
  request_in_flight_ = true;
  active_request_id_ = request_id;
  active_request_map_generation_ = request_generation;
  request_deadline_ = std::chrono::steady_clock::now() + request_timeout_;

  // Keep the callback handoff alive until both action callbacks have run. A
  // weak-only capture would let Panel destruction drop the state before a
  // delayed goal-acceptance response arrives, making it impossible to cancel
  // that accepted goal safely.
  const auto callback_state = state;
  const auto request_map_id = goal.map_id;
  // Keep the client alive until the goal response arrives. If the Panel is
  // destroyed while acceptance is pending, the callback must still be able
  // to cancel a goal that the server subsequently accepts.
  const auto action_client = action_client_;
  try {
    PlanningActionClient::SendGoalOptions options;
    options.goal_response_callback =
      [callback_state, action_client, request_id, request_generation](
        const PlanningGoalHandle::SharedPtr& goal_handle) {
        if (!goal_handle) {
          CallbackState::PendingResponse pending;
          pending.request_id = request_id;
          pending.map_generation = request_generation;
          pending.error = "The Raystar action server rejected the goal.";
          std::lock_guard<std::mutex> lock(callback_state->mutex);
          if (!callback_state->alive || callback_state->active_request_id != request_id ||
              callback_state->active_map_generation != request_generation ||
              callback_state->map_generation != request_generation) {
            return;
          }
          callback_state->pending_response = std::move(pending);
          return;
        }

        bool goal_is_current = false;
        {
          std::lock_guard<std::mutex> lock(callback_state->mutex);
          goal_is_current = callback_state->alive &&
                            callback_state->active_request_id == request_id &&
                            callback_state->active_map_generation == request_generation &&
                            callback_state->map_generation == request_generation;
          if (goal_is_current) {
            callback_state->active_planning_goal = goal_handle;
          }
        }

        // Map/topic changes, timeout, or Panel destruction may invalidate the
        // generation before the server answers the goal request. An accepted
        // stale goal must be canceled immediately instead of running unseen.
        if (!goal_is_current) {
          cancelGoalBestEffort(action_client, goal_handle);
        }
      };

    options.result_callback = [callback_state, request_id, request_generation, request_map_id](
                                const PlanningGoalHandle::WrappedResult& wrapped_result) {
      CallbackState::PendingResponse pending;
      pending.request_id = request_id;
      pending.map_generation = request_generation;
      pending.result_code = wrapped_result.code;
      if (!wrapped_result.result) {
        pending.error = "The Raystar action returned an empty result.";
      } else if (!raystar_interfaces::mapIdsEqual(wrapped_result.result->result_info.map_id,
                                                  request_map_id)) {
        pending.error = "The Raystar action result refers to a different cached map.";
      } else {
        pending.planning_response = wrapped_result.result;
      }

      std::lock_guard<std::mutex> lock(callback_state->mutex);
      if (!callback_state->alive || callback_state->active_request_id != request_id ||
          callback_state->active_map_generation != request_generation ||
          callback_state->map_generation != request_generation) {
        return;
      }
      callback_state->active_planning_goal.reset();
      callback_state->pending_response = std::move(pending);
    };

    (void)action_client->async_send_goal(goal, options);

    // A map callback may invalidate the generation while the goal request is
    // being sent. The goal-response callback cancels it if it is accepted.
    bool still_current = false;
    {
      std::lock_guard<std::mutex> lock(state->mutex);
      still_current = state->alive && state->active_request_id == request_id &&
                      state->map_generation == request_generation;
    }
    if (!still_current) {
      cancelActiveRequest();
    }
  } catch (const std::exception& exception) {
    cancelActiveRequest();
    updatePlanButtonState();
    status_label_->setText(
      boundedQString(std::string("Could not send planning request: ") + exception.what()));
  } catch (...) {
    cancelActiveRequest();
    updatePlanButtonState();
    status_label_->setText("Error: Could not send planning request.");
  }
}

void RaystarPanel::planGoalSet() {
  processCallbacks();

  const std::string selected_topic = map_topic_edit_->text().toStdString();
  if (selected_topic.empty() || selected_topic != subscribed_topic_) {
    subscribeToMap();
    return;
  }

  const std::string selected_action = goal_set_action_name_edit_->text().toStdString();
  if (selected_action.empty() || selected_action != goal_set_action_client_name_ ||
      !goal_set_action_client_) {
    configureGoalSetActionClient();
    return;
  }
  if (request_in_flight_) {
    status_label_->setText("Planning is already in progress.");
    return;
  }
  if (!goal_set_action_client_->action_server_is_ready()) {
    status_label_->setText(QString("Error: Multi-goal Action %1 is not ready")
                             .arg(boundedQString(goal_set_action_client_name_)));
    return;
  }
  if (!latest_map_ || latest_map_->data.empty()) {
    status_label_->setText("Error: No map received yet. Check map topic.");
    return;
  }
  if (latest_map_->header.frame_id.empty()) {
    status_label_->setText("Error: The received map has an empty frame_id.");
    return;
  }
  if (multi_goal_table_->rowCount() <= 0) {
    status_label_->setText("Error: Add at least one multi-goal target.");
    updatePlanButtonState();
    return;
  }

  bool start_x_ok = false;
  bool start_y_ok = false;
  const double start_x = start_x_edit_->text().toDouble(&start_x_ok);
  const double start_y = start_y_edit_->text().toDouble(&start_y_ok);
  if (!start_x_ok || !start_y_ok || !std::isfinite(start_x) || !std::isfinite(start_y)) {
    status_label_->setText("Error: Start coordinates must be finite numbers.");
    return;
  }

  GoalSetAction::Goal goal;
  goal.start.header.frame_id = latest_map_->header.frame_id;
  goal.start.pose.position.x = start_x;
  goal.start.pose.position.y = start_y;
  goal.start.pose.orientation.w = 1.0;
  goal.goals.reserve(static_cast<std::size_t>(multi_goal_table_->rowCount()));
  goal.max_path_lengths.reserve(static_cast<std::size_t>(multi_goal_table_->rowCount()));
  for (int row = 0; row < multi_goal_table_->rowCount(); ++row) {
    const auto* x_item = multi_goal_table_->item(row, 0);
    const auto* y_item = multi_goal_table_->item(row, 1);
    const auto* bound_item = multi_goal_table_->item(row, 2);
    bool x_ok = false;
    bool y_ok = false;
    bool bound_ok = false;
    const double x = x_item ? x_item->text().toDouble(&x_ok) : 0.0;
    const double y = y_item ? y_item->text().toDouble(&y_ok) : 0.0;
    const double bound = bound_item ? bound_item->text().toDouble(&bound_ok) : 0.0;
    if (!x_ok || !y_ok || !bound_ok || !std::isfinite(x) || !std::isfinite(y) ||
        !std::isfinite(bound) || bound <= 0.0) {
      status_label_->setText(
        QString("Error: Goal %1 requires finite X/Y and a finite positive budget.").arg(row));
      return;
    }
    geometry_msgs::msg::PoseStamped pose;
    pose.header.frame_id = latest_map_->header.frame_id;
    pose.pose.position.x = x;
    pose.pose.position.y = y;
    pose.pose.orientation.w = 1.0;
    goal.goals.emplace_back(std::move(pose));
    goal.max_path_lengths.push_back(bound);
  }
  goal.allow_self_crossing = allow_self_crossing_cb_->isChecked();
  goal.allow_unknown = allow_unknown_cb_->isChecked();
  goal.include_debug = request_debug_cb_->isChecked();

  auto state = callback_state_;
  if (!state) {
    status_label_->setText("Error: Panel callback state is unavailable.");
    return;
  }
  std::uint64_t request_generation = 0;
  {
    std::lock_guard<std::mutex> lock(state->mutex);
    if (!state->alive) {
      status_label_->setText("Error: Panel is shutting down.");
      return;
    }
    if (state->map_generation != map_generation_ || state->pending_map_update) {
      status_label_->setText("Map changed; wait for the new map before planning.");
      return;
    }
    request_generation = state->map_generation;
  }

  try {
    goal.map_id = raystar_interfaces::computeMapId(*latest_map_);
  } catch (const std::exception& exception) {
    status_label_->setText(
      boundedQString(std::string("Could not prepare multi-goal request: ") + exception.what()));
    return;
  } catch (...) {
    status_label_->setText("Error: Could not prepare the multi-goal request.");
    return;
  }

  std::uint64_t request_id = 0;
  {
    std::lock_guard<std::mutex> lock(state->mutex);
    if (!state->alive || state->map_generation != request_generation || state->pending_map_update) {
      status_label_->setText("Map changed; wait for the new map before planning.");
      return;
    }
    next_request_id_ = nextCounter(next_request_id_);
    request_id = next_request_id_;
    state->active_request_id = request_id;
    state->active_map_generation = request_generation;
    state->active_planning_goal.reset();
    state->active_goal_set_goal.reset();
    state->pending_response.reset();
  }

  clearResults();
  status_label_->setText(
    QString("Planning one shared tree for %1 goals...").arg(multi_goal_table_->rowCount()));
  request_in_flight_ = true;
  active_request_id_ = request_id;
  active_request_map_generation_ = request_generation;
  request_deadline_ = std::chrono::steady_clock::now() + request_timeout_;
  updatePlanButtonState();

  const auto callback_state = state;
  const auto action_client = goal_set_action_client_;
  const auto request_map_id = goal.map_id;
  const auto requested_start = goal.start;
  const auto requested_goals = goal.goals;
  const auto requested_bounds = goal.max_path_lengths;
  const bool requested_allow_self_crossing = goal.allow_self_crossing;
  const bool requested_allow_unknown = goal.allow_unknown;
  const bool requested_include_debug = goal.include_debug;
  try {
    GoalSetActionClient::SendGoalOptions options;
    options.goal_response_callback =
      [callback_state, action_client, request_id, request_generation](
        const GoalSetGoalHandle::SharedPtr& goal_handle) {
        if (!goal_handle) {
          CallbackState::PendingResponse pending;
          pending.request_id = request_id;
          pending.map_generation = request_generation;
          pending.error =
            "The Raystar multi-goal server rejected the goal (it may be busy or shutting down).";
          std::lock_guard<std::mutex> lock(callback_state->mutex);
          if (!callback_state->alive || callback_state->active_request_id != request_id ||
              callback_state->active_map_generation != request_generation ||
              callback_state->map_generation != request_generation) {
            return;
          }
          callback_state->pending_response = std::move(pending);
          return;
        }

        bool goal_is_current = false;
        {
          std::lock_guard<std::mutex> lock(callback_state->mutex);
          goal_is_current = callback_state->alive &&
                            callback_state->active_request_id == request_id &&
                            callback_state->active_map_generation == request_generation &&
                            callback_state->map_generation == request_generation;
          if (goal_is_current) {
            callback_state->active_goal_set_goal = goal_handle;
          }
        }
        if (!goal_is_current) {
          cancelGoalBestEffort(action_client, goal_handle);
        }
      };

    options.result_callback = [callback_state,
                               request_id,
                               request_generation,
                               request_map_id,
                               requested_start,
                               requested_goals,
                               requested_bounds,
                               requested_allow_self_crossing,
                               requested_allow_unknown,
                               requested_include_debug](
                                const GoalSetGoalHandle::WrappedResult& wrapped_result) {
      CallbackState::PendingResponse pending;
      pending.request_id = request_id;
      pending.map_generation = request_generation;
      pending.result_code = wrapped_result.code;
      if (!wrapped_result.result) {
        pending.error = "The Raystar multi-goal action returned an empty result.";
      } else if (!raystar_interfaces::mapIdsEqual(wrapped_result.result->result_info.map_id,
                                                  request_map_id)) {
        pending.error = "The Raystar multi-goal action result refers to a different cached map.";
      } else if (wrapped_result.result->requested_start != requested_start ||
                 wrapped_result.result->requested_allow_self_crossing !=
                   requested_allow_self_crossing ||
                 wrapped_result.result->requested_allow_unknown != requested_allow_unknown) {
        pending.error =
          "The Raystar multi-goal action result does not echo the requested start or policy.";
      } else if (wrapped_result.result->result_info.debug_requested != requested_include_debug) {
        pending.error =
          "The Raystar multi-goal action result does not echo the debug-output policy.";
      } else if (static_cast<std::size_t>(
                   wrapped_result.result->result_info.requested_goal_count) !=
                   requested_goals.size() ||
                 static_cast<std::size_t>(wrapped_result.result->result_info.returned_goal_count) !=
                   wrapped_result.result->goal_results.size()) {
        pending.error = "The Raystar multi-goal action result has inconsistent goal counts.";
      } else if ((wrapped_result.result->success ||
                  wrapped_result.result->result_info.request_satisfied) &&
                 wrapped_result.result->goal_results.size() != requested_goals.size()) {
        pending.error =
          "The Raystar multi-goal action claims completion without every goal result.";
      } else {
        bool associations_match =
          wrapped_result.result->goal_results.size() <= requested_goals.size();
        for (std::size_t i = 0;
             associations_match && i < wrapped_result.result->goal_results.size();
             ++i) {
          const auto& goal_result = wrapped_result.result->goal_results[i];
          associations_match =
            static_cast<std::size_t>(goal_result.goal_index) == i &&
            goal_result.goal == requested_goals[i] &&
            goal_result.requested_max_path_length == requested_bounds[i] &&
            goal_result.result_info.search_mode ==
              PlanningAction::Goal::SEARCH_MODE_ALL_WITHIN_LENGTH &&
            goal_result.result_info.requested_max_path_length == requested_bounds[i] &&
            goal_result.result_info.debug_requested == requested_include_debug &&
            raystar_interfaces::mapIdsEqual(goal_result.result_info.map_id, request_map_id) &&
            raystar_interfaces::environmentIdsEqual(
              goal_result.result_info.environment_id,
              wrapped_result.result->result_info.environment_id);
        }
        if (!associations_match) {
          pending.error =
            "The Raystar multi-goal action result does not preserve goal associations.";
        } else {
          pending.goal_set_response = wrapped_result.result;
        }
      }

      std::lock_guard<std::mutex> lock(callback_state->mutex);
      if (!callback_state->alive || callback_state->active_request_id != request_id ||
          callback_state->active_map_generation != request_generation ||
          callback_state->map_generation != request_generation) {
        return;
      }
      callback_state->active_goal_set_goal.reset();
      callback_state->pending_response = std::move(pending);
    };

    (void)action_client->async_send_goal(goal, options);
    bool still_current = false;
    {
      std::lock_guard<std::mutex> lock(state->mutex);
      still_current = state->alive && state->active_request_id == request_id &&
                      state->map_generation == request_generation;
    }
    if (!still_current) {
      cancelActiveRequest();
    }
  } catch (const std::exception& exception) {
    cancelActiveRequest();
    updatePlanButtonState();
    status_label_->setText(
      boundedQString(std::string("Could not send multi-goal request: ") + exception.what()));
  } catch (...) {
    cancelActiveRequest();
    updatePlanButtonState();
    status_label_->setText("Error: Could not send the multi-goal request.");
  }
}

void RaystarPanel::clearResults() {
  configurePathResultsTable(search_mode_combo_->currentIndex() == kMultiGoalWithinLengthsMode);
  results_table_->clearContents();
  results_table_->setRowCount(0);
  goal_results_table_->clearContents();
  goal_results_table_->setRowCount(0);
  node_table_->clearContents();
  node_table_->setRowCount(0);
  timing_label_->clear();
}

void RaystarPanel::cancelActiveRequest() {
  // Invalidate the logical request before touching the Action client. A
  // result callback that is already queued will then fail the generation
  // check even if cancellation races with terminal result delivery.
  PlanningGoalHandle::SharedPtr accepted_planning_goal;
  GoalSetGoalHandle::SharedPtr accepted_goal_set_goal;
  if (callback_state_) {
    std::lock_guard<std::mutex> lock(callback_state_->mutex);
    callback_state_->active_request_id = 0;
    callback_state_->active_map_generation = 0;
    callback_state_->pending_response.reset();
    accepted_planning_goal = std::move(callback_state_->active_planning_goal);
    accepted_goal_set_goal = std::move(callback_state_->active_goal_set_goal);
  }
  cancelGoalBestEffort(action_client_, accepted_planning_goal);
  cancelGoalBestEffort(goal_set_action_client_, accepted_goal_set_goal);
  request_deadline_.reset();
  request_in_flight_ = false;
  active_request_id_ = 0;
  active_request_map_generation_ = 0;
}

void RaystarPanel::processCallbacks() {
  auto state = callback_state_;
  if (!state) {
    return;
  }

  std::shared_ptr<const nav_msgs::msg::OccupancyGrid> pending_map;
  std::uint64_t pending_map_generation = 0;
  std::optional<CallbackState::PendingResponse> pending_response;
  std::deque<CallbackState::PendingPoint> pending_points;
  {
    std::lock_guard<std::mutex> lock(state->mutex);
    if (!state->alive) {
      return;
    }
    if (state->pending_map_update) {
      pending_map = std::move(state->pending_map);
      pending_map_generation = state->pending_map_generation;
      state->pending_map_generation = 0;
      state->pending_map_update = false;
    }
    if (state->pending_response) {
      pending_response = std::move(state->pending_response);
      state->pending_response.reset();
    }
    pending_points.swap(state->pending_points);
  }

  bool pending_map_is_current = false;
  if (pending_map) {
    std::lock_guard<std::mutex> lock(state->mutex);
    pending_map_is_current = state->alive && state->map_generation == pending_map_generation;
  }
  if (pending_map && pending_map_is_current && pending_map_generation != map_generation_) {
    // A new map supersedes any in-flight goal. If it was already accepted,
    // cancelActiveRequest() sends a cooperative Action cancellation request.
    if (request_in_flight_) {
      cancelActiveRequest();
    }
    latest_map_ = std::move(pending_map);
    map_generation_ = pending_map_generation;
    setCaptureMode(CaptureMode::none);
    updatePlanButtonState();
    clearResults();
    status_label_->setText(QString("Map received: %1x%2 @ %3 m/cell")
                             .arg(latest_map_->info.width)
                             .arg(latest_map_->info.height)
                             .arg(latest_map_->info.resolution, 0, 'f', 3));
  }

  for (const auto& pending : pending_points) {
    if (!latest_map_ || pending.map_generation != map_generation_) {
      continue;
    }
    if (pending.point.header.frame_id != latest_map_->header.frame_id) {
      setCaptureMode(CaptureMode::none);
      status_label_->setText(
        QString("Error: Clicked point frame '%1' does not match map frame '%2'.")
          .arg(QString::fromStdString(pending.point.header.frame_id))
          .arg(QString::fromStdString(latest_map_->header.frame_id)));
      break;
    }
    const double x = pending.point.point.x;
    const double y = pending.point.point.y;
    if (!std::isfinite(x) || !std::isfinite(y)) {
      setCaptureMode(CaptureMode::none);
      status_label_->setText("Error: Clicked point coordinates must be finite.");
      break;
    }
    if (pending.mode == CaptureMode::start_once) {
      start_x_edit_->setText(QString::number(x, 'g', 17));
      start_y_edit_->setText(QString::number(y, 'g', 17));
      setCaptureMode(CaptureMode::none);
      status_label_->setText("Start point captured.");
    } else if (pending.mode == CaptureMode::goal_once) {
      goal_x_edit_->setText(QString::number(x, 'g', 17));
      goal_y_edit_->setText(QString::number(y, 'g', 17));
      setCaptureMode(CaptureMode::none);
      status_label_->setText("Goal point captured.");
    } else if (pending.mode == CaptureMode::append_goals) {
      const int previous_goal_count = multi_goal_table_->rowCount();
      addGoalRow(x, y, max_path_length_spinbox_->value());
      if (multi_goal_table_->rowCount() > previous_goal_count) {
        status_label_->setText(
          QString("Captured goal %1; click again or toggle Capture goals to stop.")
            .arg(multi_goal_table_->rowCount() - 1));
      }
    }
  }

  if (pending_response) {
    const auto& result = *pending_response;
    bool response_is_current = false;
    {
      // Atomically claim the response against the callback-side generation.
      // A map callback that wins this lock invalidates the response before it
      // can be rendered; one that arrives later is ordered after the response
      // and its next GUI tick will clear the old result.
      std::lock_guard<std::mutex> lock(state->mutex);
      response_is_current =
        state->alive && request_in_flight_ && result.request_id == active_request_id_ &&
        result.map_generation == active_request_map_generation_ &&
        result.map_generation == map_generation_ && state->active_request_id == result.request_id &&
        state->active_map_generation == result.map_generation &&
        state->map_generation == result.map_generation;
      if (response_is_current) {
        state->active_request_id = 0;
        state->active_map_generation = 0;
      }
    }
    if (response_is_current) {
      request_in_flight_ = false;
      active_request_id_ = 0;
      active_request_map_generation_ = 0;
      request_deadline_.reset();
      updatePlanButtonState();
      try {
        displayResponse(result);
      } catch (const std::exception& exception) {
        clearResults();
        status_label_->setText(
          boundedQString(std::string("Could not display action result: ") + exception.what()));
      } catch (...) {
        clearResults();
        status_label_->setText("Error: Could not display action result.");
      }
    }
  }

  if (request_in_flight_ && request_deadline_ &&
      std::chrono::steady_clock::now() >= *request_deadline_) {
    cancelActiveRequest();
    updatePlanButtonState();
    status_label_->setText("Error: Planning goal timed out and was canceled.");
  }
}

void RaystarPanel::displayResponse(const CallbackState::PendingResponse& result) {
  if (result.planning_response) {
    displayPlanningResponse(result);
    return;
  }
  if (result.goal_set_response) {
    displayGoalSetResponse(result);
    return;
  }
  clearResults();
  status_label_->setText(QStringLiteral("Failed: ") + boundedQString(result.error));
}

void RaystarPanel::displayPlanningResponse(const CallbackState::PendingResponse& result) {
  if (!result.planning_response) {
    clearResults();
    status_label_->setText(QStringLiteral("Failed: ") + boundedQString(result.error));
    return;
  }

  configurePathResultsTable(false);
  goal_results_group_->setVisible(false);
  const auto& response = *result.planning_response;
  const std::size_t path_count = response.path_results.size();
  const int displayed_path_count = boundedRowCount(path_count);
  results_table_->setRowCount(displayed_path_count);
  for (int i = 0; i < displayed_path_count; ++i) {
    const auto index = static_cast<std::size_t>(i);
    const auto& path_result = response.path_results[index];
    results_table_->setItem(
      i, 0, new QTableWidgetItem(QString::number(static_cast<qlonglong>(index + 1))));
    results_table_->setItem(i, 1, new QTableWidgetItem(QString::number(path_result.cost, 'f', 2)));
    results_table_->setItem(i,
                            2,
                            new QTableWidgetItem(QString::number(
                              static_cast<qulonglong>(path_result.path.poses.size()))));
  }
  results_table_->resizeColumnsToContents();

  const std::size_t debug_count = response.debug_nodes.size();
  const int displayed_debug_count = boundedRowCount(debug_count);
  node_table_->setRowCount(displayed_debug_count);
  for (int i = 0; i < displayed_debug_count; ++i) {
    const auto index = static_cast<std::size_t>(i);
    const auto& debug_node = response.debug_nodes[index];
    node_table_->setItem(
      i, 0, new QTableWidgetItem(QString::number(static_cast<qlonglong>(debug_node.index))));
    node_table_->setItem(i, 1, new QTableWidgetItem(QString::number(debug_node.g_cost, 'f', 2)));
    node_table_->setItem(i, 2, new QTableWidgetItem(QString::number(debug_node.f_cost, 'f', 2)));
  }
  node_table_->resizeColumnsToContents();

  const auto& info = response.result_info;
  const QString transport_state = transportStateText(result.result_code);
  const bool bounded_mode = info.search_mode == PlanningAction::Goal::SEARCH_MODE_ALL_WITHIN_LENGTH;
  QString status;
  if (bounded_mode) {
    status = QString("%1 / %2: returned %3 path(s) within %4 m (%5 found by Core)")
               .arg(transport_state)
               .arg(planningStatusText(info.status))
               .arg(static_cast<qulonglong>(info.returned_path_count))
               .arg(info.requested_max_path_length, 0, 'g', 10)
               .arg(static_cast<qulonglong>(info.found_path_count));
  } else {
    status = QString(
               "%1 / %2: returned %3 of %4 requested "
               "(%5 found by Core)")
               .arg(transport_state)
               .arg(planningStatusText(info.status))
               .arg(static_cast<qulonglong>(info.returned_path_count))
               .arg(static_cast<qulonglong>(info.requested_path_count))
               .arg(static_cast<qulonglong>(info.found_path_count));
  }
  status += QString(
              "\nSearch complete: %1; path output complete: %2; "
              "bound exhausted: %3; limits: %4")
              .arg(info.search_complete ? QStringLiteral("yes") : QStringLiteral("no"))
              .arg(info.output_complete ? QStringLiteral("yes") : QStringLiteral("no"))
              .arg(info.cost_bound_exhausted ? QStringLiteral("yes") : QStringLiteral("no"))
              .arg(planningLimitsText(info.limits_reached));
  if (!response.message.empty()) {
    status += QStringLiteral("\n") + boundedQString(response.message);
  }
  if (info.returned_path_count != path_count) {
    status += QStringLiteral(
      "\nWarning: returned_path_count does not match "
      "path_results.size().");
  }
  if (response.success != !response.path_results.empty()) {
    status += QStringLiteral(
      "\nWarning: success does not match the structured "
      "path result count.");
  }
  if (info.found_path_count < info.returned_path_count) {
    status += QStringLiteral("\nWarning: returned path count exceeds Core found count.");
  }
  if (path_count > kMaxUiRows || debug_count > kMaxUiRows) {
    status += QStringLiteral("\nWarning: the table view was truncated for safety.");
  }
  status_label_->setText(status);
  const QString debug_summary =
    !info.debug_requested
      ? QStringLiteral("not requested")
      : QString("%1%2")
          .arg(static_cast<qulonglong>(debug_count))
          .arg(info.debug_output_complete ? QString() : QStringLiteral(" (truncated)"));
  timing_label_->setText(QString("Map: %1 ms | Search: %2 ms | Expanded: %3 | Debug: %4")
                           .arg(info.map_time_ms, 0, 'f', 1)
                           .arg(info.plan_time_ms, 0, 'f', 1)
                           .arg(static_cast<qulonglong>(info.expanded_nodes))
                           .arg(debug_summary));
}

void RaystarPanel::displayGoalSetResponse(const CallbackState::PendingResponse& result) {
  if (!result.goal_set_response) {
    clearResults();
    status_label_->setText(QStringLiteral("Failed: ") + boundedQString(result.error));
    return;
  }

  configurePathResultsTable(true);
  goal_results_group_->setVisible(true);
  const auto& response = *result.goal_set_response;
  const auto& aggregate = response.result_info;
  const std::size_t goal_count = response.goal_results.size();
  const int displayed_goal_count = boundedRowCount(goal_count);
  goal_results_table_->setRowCount(displayed_goal_count);

  std::size_t path_count = 0;
  bool per_goal_counts_consistent = true;
  bool per_goal_certificates_consistent = true;
  bool path_contracts_consistent = true;
  for (int row = 0; row < displayed_goal_count; ++row) {
    const auto& goal_result = response.goal_results[static_cast<std::size_t>(row)];
    const auto& info = goal_result.result_info;
    const bool complete = info.cost_bound_exhausted && info.output_complete;
    goal_results_table_->setItem(
      row,
      0,
      new QTableWidgetItem(QString::number(static_cast<qulonglong>(goal_result.goal_index))));
    goal_results_table_->setItem(
      row, 1, new QTableWidgetItem(QString::number(goal_result.goal.pose.position.x, 'g', 10)));
    goal_results_table_->setItem(
      row, 2, new QTableWidgetItem(QString::number(goal_result.goal.pose.position.y, 'g', 10)));
    goal_results_table_->setItem(
      row,
      3,
      new QTableWidgetItem(QString::number(goal_result.requested_max_path_length, 'g', 10)));
    QString goal_status = planningStatusText(info.status);
    if (info.limits_reached != PlanningResultInfo::LIMIT_NONE) {
      goal_status += QString(" (%1)").arg(planningLimitsText(info.limits_reached));
    }
    goal_results_table_->setItem(row, 4, new QTableWidgetItem(goal_status));
    goal_results_table_->setItem(
      row, 5, new QTableWidgetItem(complete ? QStringLiteral("yes") : QStringLiteral("no")));
    goal_results_table_->setItem(
      row,
      6,
      new QTableWidgetItem(QString("%1 / %2")
                             .arg(static_cast<qulonglong>(info.returned_path_count))
                             .arg(static_cast<qulonglong>(info.found_path_count))));

    per_goal_counts_consistent =
      per_goal_counts_consistent &&
      static_cast<std::size_t>(info.returned_path_count) == goal_result.path_results.size() &&
      info.found_path_count >= info.returned_path_count &&
      goal_result.success == !goal_result.path_results.empty();
    per_goal_certificates_consistent =
      per_goal_certificates_consistent && info.request_satisfied == complete;
    double previous_cost = -std::numeric_limits<double>::infinity();
    for (const auto& path_result : goal_result.path_results) {
      path_contracts_consistent = path_contracts_consistent && std::isfinite(path_result.cost) &&
                                  path_result.cost <= goal_result.requested_max_path_length &&
                                  path_result.cost >= previous_cost;
      previous_cost = path_result.cost;
    }
    if (path_count > std::numeric_limits<std::size_t>::max() - goal_result.path_results.size()) {
      path_count = std::numeric_limits<std::size_t>::max();
    } else {
      path_count += goal_result.path_results.size();
    }
  }
  goal_results_table_->resizeColumnsToContents();

  const int displayed_path_count = boundedRowCount(path_count);
  results_table_->setRowCount(displayed_path_count);
  int output_row = 0;
  std::size_t marker_index = 0;
  for (const auto& goal_result : response.goal_results) {
    for (std::size_t path_index = 0; path_index < goal_result.path_results.size(); ++path_index) {
      ++marker_index;
      if (output_row >= displayed_path_count) {
        continue;
      }
      const auto& path_result = goal_result.path_results[path_index];
      results_table_->setItem(
        output_row,
        0,
        new QTableWidgetItem(QString::number(static_cast<qulonglong>(goal_result.goal_index))));
      results_table_->setItem(
        output_row,
        1,
        new QTableWidgetItem(QString::number(static_cast<qulonglong>(path_index + 1))));
      results_table_->setItem(
        output_row,
        2,
        new QTableWidgetItem(response.success
                               ? QString("path_%1").arg(static_cast<qulonglong>(marker_index))
                               : QStringLiteral("not published")));
      results_table_->setItem(
        output_row, 3, new QTableWidgetItem(QString::number(path_result.cost, 'g', 10)));
      results_table_->setItem(output_row,
                              4,
                              new QTableWidgetItem(QString::number(
                                static_cast<qulonglong>(path_result.path.poses.size()))));
      ++output_row;
    }
  }
  results_table_->resizeColumnsToContents();

  const std::size_t debug_count = response.debug_nodes.size();
  const int displayed_debug_count = boundedRowCount(debug_count);
  node_table_->setRowCount(displayed_debug_count);
  for (int row = 0; row < displayed_debug_count; ++row) {
    const auto& debug_node = response.debug_nodes[static_cast<std::size_t>(row)];
    node_table_->setItem(
      row, 0, new QTableWidgetItem(QString::number(static_cast<qlonglong>(debug_node.index))));
    node_table_->setItem(row, 1, new QTableWidgetItem(QString::number(debug_node.g_cost, 'f', 2)));
    node_table_->setItem(row, 2, new QTableWidgetItem(QString::number(debug_node.f_cost, 'f', 2)));
  }
  node_table_->resizeColumnsToContents();

  QString status = QString(
                     "%1 / %2: completed %3 of %4 goals; %5 goal(s) have paths; returned %6 "
                     "path(s) (%7 found by Core)")
                     .arg(transportStateText(result.result_code))
                     .arg(planningStatusText(aggregate.status))
                     .arg(static_cast<qulonglong>(aggregate.completed_goal_count))
                     .arg(static_cast<qulonglong>(aggregate.requested_goal_count))
                     .arg(static_cast<qulonglong>(aggregate.goals_with_paths))
                     .arg(static_cast<qulonglong>(aggregate.returned_path_count))
                     .arg(static_cast<qulonglong>(aggregate.found_path_count));
  status +=
    QString("\nRequest satisfied: %1; search complete: %2; path output complete: %3; limits: %4")
      .arg(aggregate.request_satisfied ? QStringLiteral("yes") : QStringLiteral("no"))
      .arg(aggregate.search_complete ? QStringLiteral("yes") : QStringLiteral("no"))
      .arg(aggregate.output_complete ? QStringLiteral("yes") : QStringLiteral("no"))
      .arg(planningLimitsText(aggregate.limits_reached));
  if (!response.message.empty()) {
    status += QStringLiteral("\n") + boundedQString(response.message);
  }
  if (!response.success && path_count > 0) {
    status +=
      QStringLiteral("\nVisualization markers are not published for partial multi-goal results.");
  }
  if (response.success != aggregate.request_satisfied) {
    status += QStringLiteral("\nWarning: aggregate success does not match request_satisfied.");
  }
  if (static_cast<std::size_t>(aggregate.returned_goal_count) != goal_count ||
      static_cast<std::size_t>(aggregate.returned_path_count) != path_count) {
    status += QStringLiteral("\nWarning: aggregate returned counts do not match payload sizes.");
  }
  if (aggregate.found_path_count < aggregate.returned_path_count || !per_goal_counts_consistent) {
    status += QStringLiteral("\nWarning: one or more per-goal path counts are inconsistent.");
  }
  if (!per_goal_certificates_consistent) {
    status += QStringLiteral("\nWarning: one or more per-goal completion certificates disagree.");
  }
  if (!path_contracts_consistent) {
    status += QStringLiteral(
      "\nWarning: a path cost is non-finite, out of budget, or not ordered within its goal.");
  }
  if (goal_count > kMaxUiRows || path_count > kMaxUiRows || debug_count > kMaxUiRows) {
    status += QStringLiteral("\nWarning: the table view was truncated for safety.");
  }
  status_label_->setText(status);

  const QString debug_summary =
    !aggregate.debug_requested
      ? QStringLiteral("not requested")
      : QString("%1%2")
          .arg(static_cast<qulonglong>(debug_count))
          .arg(aggregate.debug_output_complete ? QString() : QStringLiteral(" (truncated)"));
  timing_label_->setText(
    QString("Map: %1 ms | Shared-tree search: %2 ms | Expanded: %3 | Debug: %4")
      .arg(aggregate.map_time_ms, 0, 'f', 1)
      .arg(aggregate.plan_time_ms, 0, 'f', 1)
      .arg(static_cast<qulonglong>(aggregate.expanded_nodes))
      .arg(debug_summary));
}

void RaystarPanel::save(rviz_common::Config config) const {
  rviz_common::Panel::save(config);
  config.mapSetValue("action_name", action_name_edit_->text());
  config.mapSetValue("goal_set_action_name", goal_set_action_name_edit_->text());
  config.mapSetValue("map_topic", map_topic_edit_->text());
  config.mapSetValue("clicked_point_topic", clicked_point_topic_edit_->text());
  config.mapSetValue("start_x", start_x_edit_->text());
  config.mapSetValue("start_y", start_y_edit_->text());
  config.mapSetValue("goal_x", goal_x_edit_->text());
  config.mapSetValue("goal_y", goal_y_edit_->text());
  config.mapSetValue("k", k_spinbox_->value());
  config.mapSetValue("search_mode", search_mode_combo_->currentIndex());
  config.mapSetValue("max_path_length",
                     QString::number(max_path_length_spinbox_->value(), 'g', 17));
  config.mapSetValue("allow_self_crossing", allow_self_crossing_cb_->isChecked());
  config.mapSetValue("allow_unknown", allow_unknown_cb_->isChecked());
  config.mapSetValue("request_debug", request_debug_cb_->isChecked());
  auto goals_config = config.mapMakeChild("multi_goals");
  for (int row = 0; row < multi_goal_table_->rowCount(); ++row) {
    auto goal_config = goals_config.listAppendNew();
    const auto* x_item = multi_goal_table_->item(row, 0);
    const auto* y_item = multi_goal_table_->item(row, 1);
    const auto* bound_item = multi_goal_table_->item(row, 2);
    goal_config.mapSetValue("x", x_item ? x_item->text() : QStringLiteral("0"));
    goal_config.mapSetValue("y", y_item ? y_item->text() : QStringLiteral("0"));
    goal_config.mapSetValue("max_path_length",
                            bound_item ? bound_item->text() : QStringLiteral("0"));
  }
}

void RaystarPanel::load(const rviz_common::Config& config) {
  rviz_common::Panel::load(config);
  bool action_changed = false;
  bool goal_set_action_changed = false;
  bool topic_changed = false;
  bool clicked_point_topic_changed = false;
  QString val;
  if (config.mapGetString("action_name", &val)) {
    action_changed = (action_name_edit_->text() != val);
    action_name_edit_->setText(val);
  }
  if (config.mapGetString("goal_set_action_name", &val)) {
    goal_set_action_changed = (goal_set_action_name_edit_->text() != val);
    goal_set_action_name_edit_->setText(val);
  }
  if (config.mapGetString("map_topic", &val)) {
    topic_changed = (map_topic_edit_->text() != val);
    map_topic_edit_->setText(val);
  }
  if (config.mapGetString("clicked_point_topic", &val)) {
    clicked_point_topic_changed = (clicked_point_topic_edit_->text() != val);
    clicked_point_topic_edit_->setText(val);
  }
  if (config.mapGetString("start_x", &val))
    start_x_edit_->setText(val);
  if (config.mapGetString("start_y", &val))
    start_y_edit_->setText(val);
  if (config.mapGetString("goal_x", &val))
    goal_x_edit_->setText(val);
  if (config.mapGetString("goal_y", &val))
    goal_y_edit_->setText(val);
  int k_val;
  if (config.mapGetInt("k", &k_val))
    k_spinbox_->setValue(k_val);
  int search_mode = 0;
  if (config.mapGetInt("search_mode", &search_mode))
    search_mode_combo_->setCurrentIndex(
      search_mode >= kTopKMode && search_mode <= kMultiGoalWithinLengthsMode ? search_mode
                                                                             : kTopKMode);
  if (config.mapGetString("max_path_length", &val)) {
    bool length_ok = false;
    const double length = val.toDouble(&length_ok);
    if (length_ok && std::isfinite(length) && length > 0.0)
      max_path_length_spinbox_->setValue(length);
  }
  bool bool_val = false;
  if (config.mapGetBool("allow_self_crossing", &bool_val)) {
    allow_self_crossing_cb_->setChecked(bool_val);
  }
  if (config.mapGetBool("allow_unknown", &bool_val)) {
    allow_unknown_cb_->setChecked(bool_val);
  }
  if (config.mapGetBool("request_debug", &bool_val)) {
    request_debug_cb_->setChecked(bool_val);
  }
  const auto goals_config = config.mapGetChild("multi_goals");
  const bool has_multi_goal_config = goals_config.isValid();
  if (has_multi_goal_config) {
    multi_goal_table_->setRowCount(0);
    const int goal_count = std::min(goals_config.listLength(), static_cast<int>(kMaxUiGoals));
    for (int row = 0; row < goal_count; ++row) {
      const auto goal_config = goals_config.listChildAt(row);
      QString x_text;
      QString y_text;
      QString bound_text;
      bool x_ok = false;
      bool y_ok = false;
      bool bound_ok = false;
      if (!goal_config.mapGetString("x", &x_text) || !goal_config.mapGetString("y", &y_text) ||
          !goal_config.mapGetString("max_path_length", &bound_text)) {
        continue;
      }
      const double x = x_text.toDouble(&x_ok);
      const double y = y_text.toDouble(&y_ok);
      const double bound = bound_text.toDouble(&bound_ok);
      if (!x_ok || !y_ok || !bound_ok || !std::isfinite(x) || !std::isfinite(y) ||
          !std::isfinite(bound) || bound <= 0.0) {
        continue;
      }
      addGoalRow(x, y, bound);
    }
  } else if (multi_goal_table_->rowCount() > 0) {
    // A pre-multi-goal RViz configuration naturally seeds the first row from
    // its legacy single-goal coordinates and bounded-search value.
    multi_goal_table_->item(0, 0)->setText(goal_x_edit_->text());
    multi_goal_table_->item(0, 1)->setText(goal_y_edit_->text());
    multi_goal_table_->item(0, 2)->setText(
      QString::number(max_path_length_spinbox_->value(), 'g', 17));
  }
  if (action_changed && ros_initialized_) {
    configureActionClient();
  }
  if (goal_set_action_changed && ros_initialized_) {
    configureGoalSetActionClient();
  }
  if (topic_changed && ros_initialized_) {
    subscribeToMap();
  }
  if (clicked_point_topic_changed && ros_initialized_) {
    subscribeToClickedPoint();
  }
  onSearchModeChanged(search_mode_combo_->currentIndex());
  updatePlanButtonState();
}

}  // namespace raystar_rviz_plugins
