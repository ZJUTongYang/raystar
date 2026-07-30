#include "raystar_rviz_plugins/raystar_panel.h"
#include <pluginlib/class_list_macros.hpp>
#include <rviz_common/display_context.hpp>
#include <QHeaderView>
#include <QStringList>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <string>
#include <utility>

PLUGINLIB_EXPORT_CLASS(raystar_rviz_plugins::RaystarPanel, rviz_common::Panel)

namespace raystar_rviz_plugins {

namespace {

constexpr std::size_t kMaxUiRows = 10000;
constexpr std::size_t kMaxUiMessageBytes = 2048;
constexpr auto kDefaultRequestTimeout = std::chrono::seconds(60);
constexpr char kDefaultPlanningActionName[] = "/raystar/plan_paths";

using PlanningAction = raystar_interfaces::action::PlanRaystarPaths;
using PlanningGoalHandle = rclcpp_action::ClientGoalHandle<PlanningAction>;
using PlanningActionClient = rclcpp_action::Client<PlanningAction>;
using PlanningResultInfo = raystar_interfaces::msg::PlanningResultInfo;

void cancelGoalBestEffort(const PlanningActionClient::SharedPtr& client,
                          const PlanningGoalHandle::SharedPtr& goal_handle) noexcept {
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
  , allow_self_crossing_cb_(new QCheckBox("Allow Self Crossing"))
  , allow_unknown_cb_(new QCheckBox("Allow Unknown"))
  , request_debug_cb_(new QCheckBox("Request Debug"))
  , action_name_edit_(new QLineEdit(kDefaultPlanningActionName))
  , map_topic_edit_(new QLineEdit("/map"))
  , plan_button_(new QPushButton("Plan"))
  , status_label_(new QLabel("Ready"))
  , results_table_(new QTableWidget())
  , node_table_(new QTableWidget())
  , timing_label_(new QLabel(""))
  , callback_state_(std::make_shared<CallbackState>())
  , callback_timer_(new QTimer(this))
  , request_timeout_(request_timeout > std::chrono::milliseconds::zero() ? request_timeout
                                                                         : kDefaultRequestTimeout) {
  k_spinbox_->setRange(1, 100);
  k_spinbox_->setValue(5);
  start_x_edit_->setObjectName("start_x_edit");
  start_y_edit_->setObjectName("start_y_edit");
  goal_x_edit_->setObjectName("goal_x_edit");
  goal_y_edit_->setObjectName("goal_y_edit");
  k_spinbox_->setObjectName("k_spinbox");
  allow_self_crossing_cb_->setObjectName("allow_self_crossing_checkbox");
  allow_unknown_cb_->setObjectName("allow_unknown_checkbox");
  request_debug_cb_->setObjectName("request_debug_checkbox");
  action_name_edit_->setObjectName("action_name_edit");
  map_topic_edit_->setObjectName("map_topic_edit");
  results_table_->setColumnCount(3);
  results_table_->setHorizontalHeaderLabels({"Path", "Cost (m)", "Waypoints"});
  results_table_->horizontalHeader()->setStretchLastSection(true);
  node_table_->setColumnCount(3);
  node_table_->setHorizontalHeaderLabels({"Node", "G-cost (m)", "F-cost (m)"});
  node_table_->horizontalHeader()->setStretchLastSection(true);
  plan_button_->setEnabled(false);
  plan_button_->setObjectName("plan_button");
  status_label_->setObjectName("status_label");
  results_table_->setObjectName("results_table");
  node_table_->setObjectName("node_table");
  timing_label_->setObjectName("timing_label");
  status_label_->setWordWrap(true);
  timing_label_->setWordWrap(true);
  setupUi();

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
    callback_state_->active_request_id = 0;
  }
  cancelActiveRequest();
  map_sub_.reset();
  action_client_.reset();
  callback_state_.reset();
}

void RaystarPanel::setupUi() {
  auto* main_layout = new QVBoxLayout;

  auto* action_group = new QGroupBox("Planner Action");
  auto* action_layout = new QHBoxLayout;
  action_layout->addWidget(new QLabel("Name:"));
  action_layout->addWidget(action_name_edit_);
  action_group->setLayout(action_layout);
  main_layout->addWidget(action_group);

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
  param_layout->addWidget(request_debug_cb_);
  main_layout->addLayout(param_layout);

  main_layout->addWidget(plan_button_);
  main_layout->addWidget(status_label_);
  main_layout->addWidget(results_table_);
  main_layout->addWidget(new QLabel("Tree nodes:"));
  main_layout->addWidget(node_table_);
  main_layout->addWidget(timing_label_);
  main_layout->addStretch();

  setLayout(main_layout);

  connect(plan_button_, &QPushButton::clicked, this, &RaystarPanel::onPlanClicked);
  connect(map_topic_edit_, &QLineEdit::textChanged, this, &RaystarPanel::onMapTopicChanged);
  connect(map_topic_edit_, &QLineEdit::editingFinished, this, &RaystarPanel::subscribeToMap);
  connect(action_name_edit_, &QLineEdit::textChanged, this, &RaystarPanel::onActionNameChanged);
  connect(
    action_name_edit_, &QLineEdit::editingFinished, this, &RaystarPanel::configureActionClient);
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
  configureActionClient();
}

void RaystarPanel::updatePlanButtonState() {
  const bool action_is_current =
    action_client_ && action_name_edit_->text().toStdString() == action_client_name_;
  const bool map_is_current = latest_map_ && !latest_map_->data.empty() &&
                              map_topic_edit_->text().toStdString() == subscribed_topic_;
  plan_button_->setEnabled(ros_initialized_ && action_is_current && map_is_current &&
                           !request_in_flight_);
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

std::uint64_t RaystarPanel::invalidateMapState() {
  map_generation_ = nextCounter(map_generation_);
  cancelActiveRequest();
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
      });
    subscribed_topic_ = topic;
    status_label_->setText(QString("Waiting for map on %1").arg(boundedQString(topic)));
  } catch (const std::exception&) {
    map_sub_.reset();
    status_label_->setText("Error: Could not subscribe to the selected map topic.");
  }
}

void RaystarPanel::onPlanClicked() {
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
    goal.k = k_spinbox_->value();
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
    state->active_goal.reset();
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
            callback_state->active_goal = goal_handle;
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
        pending.response = wrapped_result.result;
      }

      std::lock_guard<std::mutex> lock(callback_state->mutex);
      if (!callback_state->alive || callback_state->active_request_id != request_id ||
          callback_state->active_map_generation != request_generation ||
          callback_state->map_generation != request_generation) {
        return;
      }
      callback_state->active_goal.reset();
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

void RaystarPanel::clearResults() {
  results_table_->clearContents();
  results_table_->setRowCount(0);
  node_table_->clearContents();
  node_table_->setRowCount(0);
  timing_label_->clear();
}

void RaystarPanel::cancelActiveRequest() {
  // Invalidate the logical request before touching the Action client. A
  // result callback that is already queued will then fail the generation
  // check even if cancellation races with terminal result delivery.
  PlanningGoalHandle::SharedPtr accepted_goal;
  if (callback_state_) {
    std::lock_guard<std::mutex> lock(callback_state_->mutex);
    callback_state_->active_request_id = 0;
    callback_state_->active_map_generation = 0;
    callback_state_->pending_response.reset();
    accepted_goal = std::move(callback_state_->active_goal);
  }
  cancelGoalBestEffort(action_client_, accepted_goal);
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
    updatePlanButtonState();
    clearResults();
    status_label_->setText(QString("Map received: %1x%2 @ %3 m/cell")
                             .arg(latest_map_->info.width)
                             .arg(latest_map_->info.height)
                             .arg(latest_map_->info.resolution, 0, 'f', 3));
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
  if (!result.response) {
    clearResults();
    status_label_->setText(QStringLiteral("Failed: ") + boundedQString(result.error));
    return;
  }

  const auto& response = *result.response;
  const bool action_succeeded = result.result_code == rclcpp_action::ResultCode::SUCCEEDED;
  const bool action_canceled = result.result_code == rclcpp_action::ResultCode::CANCELED;
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
  const QString transport_state =
    action_succeeded ? QStringLiteral("SUCCEEDED")
                     : (action_canceled ? QStringLiteral("CANCELED")
                                        : (result.result_code == rclcpp_action::ResultCode::ABORTED
                                             ? QStringLiteral("ABORTED")
                                             : QStringLiteral("UNKNOWN")));
  QString status = QString(
                     "%1 / %2: returned %3 of %4 requested "
                     "(%5 found by Core)")
                     .arg(transport_state)
                     .arg(planningStatusText(info.status))
                     .arg(static_cast<qulonglong>(info.returned_path_count))
                     .arg(static_cast<qulonglong>(info.requested_path_count))
                     .arg(static_cast<qulonglong>(info.found_path_count));
  status += QString(
              "\nSearch complete: %1; path output complete: %2; "
              "limits: %3")
              .arg(info.search_complete ? QStringLiteral("yes") : QStringLiteral("no"))
              .arg(info.output_complete ? QStringLiteral("yes") : QStringLiteral("no"))
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

void RaystarPanel::save(rviz_common::Config config) const {
  rviz_common::Panel::save(config);
  config.mapSetValue("action_name", action_name_edit_->text());
  config.mapSetValue("map_topic", map_topic_edit_->text());
  config.mapSetValue("start_x", start_x_edit_->text());
  config.mapSetValue("start_y", start_y_edit_->text());
  config.mapSetValue("goal_x", goal_x_edit_->text());
  config.mapSetValue("goal_y", goal_y_edit_->text());
  config.mapSetValue("k", k_spinbox_->value());
  config.mapSetValue("allow_self_crossing", allow_self_crossing_cb_->isChecked());
  config.mapSetValue("allow_unknown", allow_unknown_cb_->isChecked());
  config.mapSetValue("request_debug", request_debug_cb_->isChecked());
}

void RaystarPanel::load(const rviz_common::Config& config) {
  rviz_common::Panel::load(config);
  bool action_changed = false;
  bool topic_changed = false;
  QString val;
  if (config.mapGetString("action_name", &val)) {
    action_changed = (action_name_edit_->text() != val);
    action_name_edit_->setText(val);
  }
  if (config.mapGetString("map_topic", &val)) {
    topic_changed = (map_topic_edit_->text() != val);
    map_topic_edit_->setText(val);
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
  if (action_changed && ros_initialized_) {
    configureActionClient();
  }
  if (topic_changed && ros_initialized_) {
    subscribeToMap();
  }
  updatePlanButtonState();
}

}  // namespace raystar_rviz_plugins
