#include "raystar_node.h"
#include "conservative_path_length.h"
#include "metric_bound_search.h"
#include "published_path_order.h"

#include <std_msgs/msg/color_rgba.hpp>
#include <raystar_interfaces/msg/debug_node.hpp>
#include <raystar_interfaces/msg/path_result.hpp>
#include <raystar_interfaces/msg/planning_result_info.hpp>
#include <rcl_interfaces/msg/integer_range.hpp>
#include <rcl_interfaces/msg/parameter_descriptor.hpp>
#include <rcl_interfaces/msg/set_parameters_result.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <iomanip>
#include <type_traits>
#include <utility>
#include <exception>

#include "raystar_node_detail.h"

namespace raystar {

using namespace node_impl;

void RaystarNode::handleActionAccepted(const std::shared_ptr<PlanGoalHandle> goal_handle) {
  std::shared_ptr<std::atomic<bool>> cancel_requested;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (goal_handle && active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
      cancel_requested = active_goal_cancel_;
    }
  }

  if (!goal_handle || !cancel_requested) {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    active_goal_reserved_ = false;
    active_goal_cancel_.reset();
    planning_busy_.store(false, std::memory_order_release);
    return;
  }

  bool queued = false;
  {
    std::lock_guard<std::mutex> lock(action_worker_mutex_);
    if (!stop_action_worker_ && !pending_action_job_) {
      pending_action_job_.emplace(ActionJob{goal_handle, cancel_requested});
      queued = true;
    }
  }
  if (queued) {
    action_worker_cv_.notify_one();
    return;
  }

  // Release admission before doing any best-effort result allocation. This
  // branch is only reachable during shutdown or an internal invariant
  // violation; even an allocation failure must not leave the planner busy.
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
      active_goal_reserved_ = false;
      active_goal_cancel_.reset();
    }
    planning_busy_.store(false, std::memory_order_release);
  }
  try {
    auto result = std::make_shared<PlanAction::Result>();
    const auto rejected_goal = goal_handle ? goal_handle->get_goal() : nullptr;
    if (rejected_goal) {
      initializePlanningResponse(*result,
                                 rejected_goal->search_mode,
                                 rejected_goal->k,
                                 rejected_goal->max_path_length,
                                 rejected_goal->map_id,
                                 rejected_goal->include_debug);
    } else {
      resetPlanningResponse(*result);
    }
    result->result_info.status = PlanningResultInfo::STATUS_FAILED;
    result->message = "Raystar Action worker is unavailable";
    goal_handle->abort(result);
  } catch (...) {}
}

rclcpp_action::GoalResponse RaystarNode::handleGoalSetActionGoal(
  const rclcpp_action::GoalUUID& uuid, std::shared_ptr<const GoalSetAction::Goal> goal) {
  if (!goal || shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
    return rclcpp_action::GoalResponse::REJECT;

  std::shared_ptr<std::atomic<bool>> cancel_requested;
  try {
    cancel_requested = std::make_shared<std::atomic<bool>>(false);
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Rejecting goal-set Action: could not allocate cancellation state");
    return rclcpp_action::GoalResponse::REJECT;
  }
  bool expected_idle = false;
  if (!planning_busy_.compare_exchange_strong(expected_idle, true, std::memory_order_acq_rel)) {
    RCLCPP_WARN(get_logger(), "Rejecting goal-set Action because the capacity-one planner is busy");
    return rclcpp_action::GoalResponse::REJECT;
  }
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    active_goal_id_ = uuid;
    active_goal_reserved_ = true;
    active_goal_cancel_ = std::move(cancel_requested);
  }
  return rclcpp_action::GoalResponse::ACCEPT_AND_EXECUTE;
}

rclcpp_action::CancelResponse RaystarNode::handleGoalSetActionCancel(
  const std::shared_ptr<GoalSetGoalHandle> goal_handle) {
  if (!goal_handle)
    return rclcpp_action::CancelResponse::REJECT;
  std::lock_guard<std::mutex> lock(action_state_mutex_);
  if (!active_goal_reserved_ || active_goal_id_ != goal_handle->get_goal_id() ||
      !active_goal_cancel_) {
    return rclcpp_action::CancelResponse::REJECT;
  }
  active_goal_cancel_->store(true, std::memory_order_release);
  return rclcpp_action::CancelResponse::ACCEPT;
}

void RaystarNode::handleGoalSetActionAccepted(
  const std::shared_ptr<GoalSetGoalHandle> goal_handle) {
  std::shared_ptr<std::atomic<bool>> cancel_requested;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (goal_handle && active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id())
      cancel_requested = active_goal_cancel_;
  }
  if (!goal_handle || !cancel_requested) {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    active_goal_reserved_ = false;
    active_goal_cancel_.reset();
    planning_busy_.store(false, std::memory_order_release);
    return;
  }

  bool queued = false;
  {
    std::lock_guard<std::mutex> lock(action_worker_mutex_);
    if (!stop_action_worker_ && !pending_action_job_) {
      pending_action_job_.emplace(GoalSetActionJob{goal_handle, cancel_requested});
      queued = true;
    }
  }
  if (queued) {
    action_worker_cv_.notify_one();
    return;
  }

  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
      active_goal_reserved_ = false;
      active_goal_cancel_.reset();
    }
    planning_busy_.store(false, std::memory_order_release);
  }
  try {
    auto result = std::make_shared<GoalSetAction::Result>();
    result->result_info.map_id = goal_handle->get_goal()->map_id;
    result->result_info.status = PlanningResultInfo::STATUS_FAILED;
    result->message = "Raystar Action worker is unavailable";
    goal_handle->abort(result);
  } catch (...) {}
}

rclcpp_action::GoalResponse RaystarNode::handleTransitionActionGoal(
  const rclcpp_action::GoalUUID& uuid, std::shared_ptr<const TransitionAction::Goal> goal) {
  if (!goal || shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
    return rclcpp_action::GoalResponse::REJECT;

  std::shared_ptr<std::atomic<bool>> cancel_requested;
  try {
    cancel_requested = std::make_shared<std::atomic<bool>>(false);
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Rejecting UPS Action: could not allocate cancellation state");
    return rclcpp_action::GoalResponse::REJECT;
  }
  bool expected_idle = false;
  if (!planning_busy_.compare_exchange_strong(expected_idle, true, std::memory_order_acq_rel)) {
    RCLCPP_WARN(get_logger(), "Rejecting UPS Action because the capacity-one planner is busy");
    return rclcpp_action::GoalResponse::REJECT;
  }
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    active_goal_id_ = uuid;
    active_goal_reserved_ = true;
    active_goal_cancel_ = std::move(cancel_requested);
  }
  return rclcpp_action::GoalResponse::ACCEPT_AND_EXECUTE;
}

rclcpp_action::CancelResponse RaystarNode::handleTransitionActionCancel(
  const std::shared_ptr<TransitionGoalHandle> goal_handle) {
  if (!goal_handle)
    return rclcpp_action::CancelResponse::REJECT;
  std::lock_guard<std::mutex> lock(action_state_mutex_);
  if (!active_goal_reserved_ || active_goal_id_ != goal_handle->get_goal_id() ||
      !active_goal_cancel_) {
    return rclcpp_action::CancelResponse::REJECT;
  }
  active_goal_cancel_->store(true, std::memory_order_release);
  return rclcpp_action::CancelResponse::ACCEPT;
}

void RaystarNode::handleTransitionActionAccepted(
  const std::shared_ptr<TransitionGoalHandle> goal_handle) {
  std::shared_ptr<std::atomic<bool>> cancel_requested;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (goal_handle && active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id())
      cancel_requested = active_goal_cancel_;
  }
  if (!goal_handle || !cancel_requested) {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    active_goal_reserved_ = false;
    active_goal_cancel_.reset();
    planning_busy_.store(false, std::memory_order_release);
    return;
  }

  bool queued = false;
  {
    std::lock_guard<std::mutex> lock(action_worker_mutex_);
    if (!stop_action_worker_ && !pending_action_job_) {
      pending_action_job_.emplace(TransitionActionJob{goal_handle, cancel_requested});
      queued = true;
    }
  }
  if (queued) {
    action_worker_cv_.notify_one();
    return;
  }

  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
      active_goal_reserved_ = false;
      active_goal_cancel_.reset();
    }
    planning_busy_.store(false, std::memory_order_release);
  }
  try {
    auto result = std::make_shared<TransitionAction::Result>();
    result->map_id = goal_handle->get_goal()->map_id;
    result->status = TransitionAction::Result::STATUS_FAILED;
    result->message = "Raystar Action worker is unavailable";
    goal_handle->abort(result);
  } catch (...) {}
}

void RaystarNode::actionWorkerLoop() noexcept {
  while (true) {
    std::optional<PendingActionJob> job;
    {
      std::unique_lock<std::mutex> lock(action_worker_mutex_);
      action_worker_cv_.wait(
        lock, [this]() { return stop_action_worker_ || pending_action_job_.has_value(); });
      if (stop_action_worker_ && !pending_action_job_)
        return;
      job.emplace(std::move(*pending_action_job_));
      pending_action_job_.reset();
    }
    std::visit(
      [this](auto&& typed_job) {
        using Job = std::decay_t<decltype(typed_job)>;
        if constexpr (std::is_same_v<Job, ActionJob>) {
          executeAction(typed_job.goal_handle, typed_job.cancel_requested);
        } else if constexpr (std::is_same_v<Job, GoalSetActionJob>) {
          executeGoalSetAction(typed_job.goal_handle, typed_job.cancel_requested);
        } else {
          executeTransitionAction(typed_job.goal_handle, typed_job.cancel_requested);
        }
      },
      *job);
  }
}

void RaystarNode::executeAction(
  const std::shared_ptr<PlanGoalHandle> goal_handle,
  const std::shared_ptr<std::atomic<bool>>& cancel_requested) noexcept {
  enum class TerminalState { succeeded, aborted, canceled };

  std::shared_ptr<PlanAction::Result> result;
  TerminalState terminal_state = TerminalState::aborted;
  const auto goal_is_canceling = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_canceling();
    } catch (...) {
      return false;
    }
  };
  const auto goal_is_active = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_active();
    } catch (...) {
      return false;
    }
  };

  try {
    result = std::make_shared<PlanAction::Result>();
    resetPlanningResponse(*result);
    const std::weak_ptr<PlanGoalHandle> weak_goal_handle(goal_handle);
    const StopPredicate stop_requested = [this, weak_goal_handle, cancel_requested]() {
      if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
        return true;
      const auto handle = weak_goal_handle.lock();
      if (!handle)
        return true;
      // Stop Core as soon as our cancel callback accepts the request.  The
      // worker separately waits for rclcpp_action's subsequent CANCELING
      // transition before it publishes the terminal result.
      return cancel_requested->load(std::memory_order_acquire);
    };

    const auto goal = goal_handle->get_goal();
    if (goal) {
      nav_msgs::msg::OccupancyGrid::ConstSharedPtr cached_map;
      std::string map_error;
      if (!resolveCachedMap(goal->map_id, cached_map, map_error)) {
        initializePlanningResponse(*result,
                                   goal->search_mode,
                                   goal->k,
                                   goal->max_path_length,
                                   goal->map_id,
                                   goal->include_debug);
        result->result_info.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
        result->message = std::move(map_error);
      } else {
        executePlanning(*goal, *result, *cached_map, goal->map_id, stop_requested);
      }
    } else {
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->message = "Raystar Action goal data is unavailable";
    }

    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok()) {
      markFailed(*result);
      if (result->message.empty())
        result->message = "Planning stopped because the Raystar node is shutting down";
      terminal_state = TerminalState::aborted;
    } else if (goal_is_canceling()) {
      markCancelled(*result);
      if (result->message.empty())
        result->message = "Planning was cancelled";
      terminal_state = TerminalState::canceled;
    } else if (result->success) {
      terminal_state = TerminalState::succeeded;
    } else {
      terminal_state = TerminalState::aborted;
    }
  } catch (const std::exception& exception) {
    try {
      if (!result)
        result = std::make_shared<PlanAction::Result>();
      resetPlanningResponsePreservingRequestMetadata(*result);
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      setBoundedExceptionMessage(result->message, exception.what());
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
    RCLCPP_ERROR(get_logger(), "Raystar Action worker failed: %s", exception.what());
  } catch (...) {
    try {
      if (!result)
        result = std::make_shared<PlanAction::Result>();
      resetPlanningResponsePreservingRequestMetadata(*result);
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->message = "Raystar Action worker failed with an unknown exception";
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
    RCLCPP_ERROR(get_logger(), "Raystar Action worker failed with an unknown exception");
  }

  // Linearize completion against handleActionCancel().  If the cancel
  // callback stored its per-goal flag first, retain the reservation while
  // rclcpp_action performs its internal transition after that callback
  // returns.  Otherwise clear the reservation under the same mutex, causing a
  // later cancel request to be rejected before this normal terminal result.
  bool explicit_cancel = false;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    explicit_cancel = cancel_requested->load(std::memory_order_acquire) &&
                      !shutting_down_.load(std::memory_order_acquire);
    if (!explicit_cancel) {
      if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
        terminal_state = TerminalState::aborted;
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
      planning_busy_.store(false, std::memory_order_release);
    }
  }

  // An accepted cancel must retain the capacity-one reservation through both
  // Marker/cache cleanup and the Action terminal transition.  Releasing busy
  // earlier would allow a new request to publish a fresh snapshot which this
  // older goal could then erase in its normal or abort-to-cancel fallback.
  ScopeExit explicit_cancel_release([&]() noexcept {
    if (!explicit_cancel)
      return;
    try {
      std::lock_guard<std::mutex> lock(action_state_mutex_);
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
    } catch (...) {
      // Admission still has to be released if middleware state is already
      // tearing down; the node destructor owns the remaining cleanup.
    }
    planning_busy_.store(false, std::memory_order_release);
  });

  if (explicit_cancel) {
    terminal_state = TerminalState::canceled;
    const auto transition_deadline = std::chrono::steady_clock::now() + std::chrono::seconds(1);
    while (rclcpp::ok() && goal_is_active() && !goal_is_canceling() &&
           std::chrono::steady_clock::now() < transition_deadline) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (!goal_is_canceling()) {
      terminal_state = TerminalState::aborted;
      if (result) {
        resetPlanningResponsePreservingRequestMetadata(*result);
        markFailed(*result);
        result->message =
          "Cancellation was accepted but the Action state transition did not complete";
      }
    }
  }

  // Allocation failure is the only case in which no typed Action result can
  // be published.  The reservation has still been released and destruction
  // of the goal handle will force a terminal state inside rclcpp_action.
  if (!result)
    return;

  if (terminal_state == TerminalState::canceled) {
    markCancelled(*result);
    if (result->message.empty())
      result->message = "Planning was cancelled";
    std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
    clearVisualizationsLocked();
  }

  try {
    if (terminal_state == TerminalState::canceled)
      goal_handle->canceled(result);
    else if (terminal_state == TerminalState::succeeded)
      goal_handle->succeed(result);
    else
      goal_handle->abort(result);
  } catch (const std::exception& exception) {
    // The bounded cancel-transition fallback can race with the transition at
    // its deadline.  If abort() lost that race, finish with the now-valid
    // canceled transition instead of leaving the goal non-terminal.
    if (terminal_state != TerminalState::canceled && goal_is_canceling()) {
      try {
        markCancelled(*result);
        result->message = "Planning was cancelled";
        {
          std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
          clearVisualizationsLocked();
        }
        goal_handle->canceled(result);
        return;
      } catch (...) {}
    }
    RCLCPP_ERROR(
      get_logger(), "Could not publish Raystar Action terminal result: %s", exception.what());
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Could not publish Raystar Action terminal result");
  }
}

void RaystarNode::executeGoalSetAction(
  const std::shared_ptr<GoalSetGoalHandle> goal_handle,
  const std::shared_ptr<std::atomic<bool>>& cancel_requested) noexcept {
  enum class TerminalState { succeeded, aborted, canceled };
  std::shared_ptr<GoalSetAction::Result> result;
  TerminalState terminal_state = TerminalState::aborted;
  const auto goal_is_canceling = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_canceling();
    } catch (...) {
      return false;
    }
  };
  const auto goal_is_active = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_active();
    } catch (...) {
      return false;
    }
  };
  const auto mark_cancelled = [](GoalSetAction::Result& output) {
    output.success = false;
    output.goal_results.clear();
    output.debug_nodes.clear();
    output.result_info.status = PlanningResultInfo::STATUS_CANCELLED;
    output.result_info.limits_reached = static_cast<uint16_t>(output.result_info.limits_reached |
                                                              PlanningResultInfo::LIMIT_CANCELLED);
    output.result_info.request_satisfied = false;
    output.result_info.search_complete = false;
    output.result_info.output_complete = false;
    output.result_info.debug_output_complete = false;
    output.result_info.returned_goal_count = 0;
    output.result_info.completed_goal_count = 0;
    output.result_info.goals_with_paths = 0;
    output.result_info.found_path_count = 0;
    output.result_info.returned_path_count = 0;
  };

  try {
    result = std::make_shared<GoalSetAction::Result>();
    const std::weak_ptr<GoalSetGoalHandle> weak_goal_handle(goal_handle);
    const StopPredicate stop_requested = [this, weak_goal_handle, cancel_requested]() {
      if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
        return true;
      if (!weak_goal_handle.lock())
        return true;
      return cancel_requested->load(std::memory_order_acquire);
    };
    const auto goal = goal_handle->get_goal();
    if (goal) {
      result->result_info.map_id = goal->map_id;
      result->result_info.requested_goal_count = boundedUint32(goal->goals.size());
      result->result_info.debug_requested = goal->include_debug;
      result->result_info.debug_output_complete = !goal->include_debug;
      nav_msgs::msg::OccupancyGrid::ConstSharedPtr cached_map;
      std::string map_error;
      if (!resolveCachedMap(goal->map_id, cached_map, map_error)) {
        result->result_info.status = PlanningResultInfo::STATUS_INVALID_REQUEST;
        result->message = std::move(map_error);
      } else {
        executeGoalSetPlanning(*goal, *result, *cached_map, goal->map_id, stop_requested);
      }
    } else {
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->message = "Raystar goal-set Action data is unavailable";
    }

    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok()) {
      result->success = false;
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->result_info.request_satisfied = false;
      result->result_info.search_complete = false;
      terminal_state = TerminalState::aborted;
    } else if (goal_is_canceling()) {
      mark_cancelled(*result);
      if (result->message.empty())
        result->message = "Planning was cancelled";
      terminal_state = TerminalState::canceled;
    } else if (result->success) {
      terminal_state = TerminalState::succeeded;
    }
  } catch (const std::exception& exception) {
    try {
      if (!result)
        result = std::make_shared<GoalSetAction::Result>();
      result->success = false;
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->goal_results.clear();
      result->debug_nodes.clear();
      setBoundedExceptionMessage(result->message, exception.what());
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
    RCLCPP_ERROR(get_logger(), "Raystar goal-set Action worker failed: %s", exception.what());
  } catch (...) {
    try {
      if (!result)
        result = std::make_shared<GoalSetAction::Result>();
      result->success = false;
      result->result_info.status = PlanningResultInfo::STATUS_FAILED;
      result->goal_results.clear();
      result->debug_nodes.clear();
      result->message = "Raystar goal-set Action worker failed with an unknown exception";
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
  }

  bool explicit_cancel = false;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    explicit_cancel = cancel_requested->load(std::memory_order_acquire) &&
                      !shutting_down_.load(std::memory_order_acquire);
    if (!explicit_cancel) {
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
      planning_busy_.store(false, std::memory_order_release);
    }
  }
  ScopeExit explicit_cancel_release([&]() noexcept {
    if (!explicit_cancel)
      return;
    try {
      std::lock_guard<std::mutex> lock(action_state_mutex_);
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
    } catch (...) {}
    planning_busy_.store(false, std::memory_order_release);
  });

  if (explicit_cancel) {
    terminal_state = TerminalState::canceled;
    const auto transition_deadline = std::chrono::steady_clock::now() + std::chrono::seconds(1);
    while (rclcpp::ok() && goal_is_active() && !goal_is_canceling() &&
           std::chrono::steady_clock::now() < transition_deadline) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (!goal_is_canceling()) {
      terminal_state = TerminalState::aborted;
      if (result) {
        result->success = false;
        result->result_info.status = PlanningResultInfo::STATUS_FAILED;
        result->message =
          "Cancellation was accepted but the Action state transition did not complete";
      }
    }
  }
  if (!result)
    return;
  if (terminal_state == TerminalState::canceled) {
    mark_cancelled(*result);
    if (result->message.empty())
      result->message = "Planning was cancelled";
    std::lock_guard<std::mutex> planner_lock(planner_cache_mutex_);
    clearVisualizationsLocked();
  }
  try {
    if (terminal_state == TerminalState::canceled)
      goal_handle->canceled(result);
    else if (terminal_state == TerminalState::succeeded)
      goal_handle->succeed(result);
    else
      goal_handle->abort(result);
  } catch (const std::exception& exception) {
    RCLCPP_ERROR(
      get_logger(), "Could not publish Raystar goal-set Action result: %s", exception.what());
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Could not publish Raystar goal-set Action result");
  }
}

void RaystarNode::executeTransitionAction(
  const std::shared_ptr<TransitionGoalHandle> goal_handle,
  const std::shared_ptr<std::atomic<bool>>& cancel_requested) noexcept {
  enum class TerminalState { succeeded, aborted, canceled };
  std::shared_ptr<TransitionAction::Result> result;
  TerminalState terminal_state = TerminalState::aborted;
  const auto goal_is_canceling = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_canceling();
    } catch (...) {
      return false;
    }
  };
  const auto goal_is_active = [this, &goal_handle]() noexcept {
    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
      return false;
    try {
      return goal_handle && goal_handle->is_active();
    } catch (...) {
      return false;
    }
  };
  const auto mark_cancelled = [](TransitionAction::Result& output) {
    output.success = false;
    output.status = TransitionAction::Result::STATUS_CANCELLED;
    if (output.message.empty())
      output.message = "UPS transition construction was cancelled";
  };

  try {
    result = std::make_shared<TransitionAction::Result>();
    const std::weak_ptr<TransitionGoalHandle> weak_goal_handle(goal_handle);
    const StopPredicate stop_requested = [this, weak_goal_handle, cancel_requested]() {
      if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok())
        return true;
      if (!weak_goal_handle.lock())
        return true;
      return cancel_requested->load(std::memory_order_acquire);
    };
    const TransitionProgressCallback progress_callback =
      [weak_goal_handle](
        std::uint32_t requested, std::uint32_t completed, const char* stage) noexcept {
        try {
          const auto handle = weak_goal_handle.lock();
          if (!handle || !handle->is_active() || !handle->is_executing())
            return;
          auto feedback = std::make_shared<TransitionAction::Feedback>();
          feedback->requested_transition_count = requested;
          feedback->completed_transition_count = completed;
          feedback->stage = stage == nullptr ? "" : stage;
          handle->publish_feedback(feedback);
        } catch (...) {
          // Feedback is best-effort. A canceled/expired handle, allocation
          // failure, or middleware exception must not abort UPS planning.
        }
      };
    const auto goal = goal_handle->get_goal();
    if (goal) {
      result->map_id = goal->map_id;
      result->requested_transition_count = boundedUint32(goal->transition_pairs.size());
      nav_msgs::msg::OccupancyGrid::ConstSharedPtr cached_map;
      std::string map_error;
      if (!resolveCachedMap(goal->map_id, cached_map, map_error)) {
        progress_callback(result->requested_transition_count, 0, "validating transition request");
        result->status = TransitionAction::Result::STATUS_INVALID_REQUEST;
        result->message = std::move(map_error);
      } else {
        executeTransitionPlanning(
          *goal, *result, *cached_map, goal->map_id, stop_requested, progress_callback);
      }
    } else {
      result->status = TransitionAction::Result::STATUS_FAILED;
      result->message = "Raystar UPS Action data is unavailable";
    }

    if (shutting_down_.load(std::memory_order_acquire) || !rclcpp::ok()) {
      result->success = false;
      result->status = TransitionAction::Result::STATUS_FAILED;
      terminal_state = TerminalState::aborted;
    } else if (goal_is_canceling()) {
      mark_cancelled(*result);
      terminal_state = TerminalState::canceled;
    } else if (result->success) {
      terminal_state = TerminalState::succeeded;
    }
  } catch (const std::exception& exception) {
    try {
      if (!result)
        result = std::make_shared<TransitionAction::Result>();
      result->success = false;
      result->status = TransitionAction::Result::STATUS_FAILED;
      result->transitions.clear();
      result->completed_transition_count = 0;
      setBoundedExceptionMessage(result->message, exception.what());
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
    RCLCPP_ERROR(get_logger(), "Raystar UPS Action worker failed: %s", exception.what());
  } catch (...) {
    try {
      if (!result)
        result = std::make_shared<TransitionAction::Result>();
      result->success = false;
      result->status = TransitionAction::Result::STATUS_FAILED;
      result->transitions.clear();
      result->completed_transition_count = 0;
      result->message = "Raystar UPS Action worker failed with an unknown exception";
    } catch (...) {}
    terminal_state = goal_is_canceling() ? TerminalState::canceled : TerminalState::aborted;
  }

  bool explicit_cancel = false;
  {
    std::lock_guard<std::mutex> lock(action_state_mutex_);
    explicit_cancel = cancel_requested->load(std::memory_order_acquire) &&
                      !shutting_down_.load(std::memory_order_acquire);
    if (!explicit_cancel) {
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
      planning_busy_.store(false, std::memory_order_release);
    }
  }
  ScopeExit explicit_cancel_release([&]() noexcept {
    if (!explicit_cancel)
      return;
    try {
      std::lock_guard<std::mutex> lock(action_state_mutex_);
      if (active_goal_reserved_ && active_goal_id_ == goal_handle->get_goal_id()) {
        active_goal_reserved_ = false;
        active_goal_cancel_.reset();
      }
    } catch (...) {}
    planning_busy_.store(false, std::memory_order_release);
  });

  if (explicit_cancel) {
    terminal_state = TerminalState::canceled;
    const auto transition_deadline = std::chrono::steady_clock::now() + std::chrono::seconds(1);
    while (rclcpp::ok() && goal_is_active() && !goal_is_canceling() &&
           std::chrono::steady_clock::now() < transition_deadline) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if (!goal_is_canceling()) {
      terminal_state = TerminalState::aborted;
      if (result) {
        result->success = false;
        result->status = TransitionAction::Result::STATUS_FAILED;
        result->message =
          "Cancellation was accepted but the Action state transition did not complete";
      }
    }
  }
  if (!result)
    return;
  if (terminal_state == TerminalState::canceled)
    mark_cancelled(*result);
  try {
    if (terminal_state == TerminalState::canceled)
      goal_handle->canceled(result);
    else if (terminal_state == TerminalState::succeeded)
      goal_handle->succeed(result);
    else
      goal_handle->abort(result);
  } catch (const std::exception& exception) {
    RCLCPP_ERROR(get_logger(), "Could not publish Raystar UPS Action result: %s", exception.what());
  } catch (...) {
    RCLCPP_ERROR(get_logger(), "Could not publish Raystar UPS Action result");
  }
}


}  // namespace raystar
