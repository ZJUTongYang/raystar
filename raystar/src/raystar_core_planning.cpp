#include <raystar/raystar_core.h>
#include <chrono>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <limits>
#include <exception>
#include <optional>

#include "exact_point_location.h"
#include "conservative_path_length.h"
#include "visibility_validation.h"

#include "raystar_core_detail.h"

namespace raystar {

using namespace core_impl;

HomotopyShorteningResult RaystarCore::shortenWithinHomotopy(const Polymap& polymap,
                                                            const PathSolution& alpha_a,
                                                            const PathSolution& alpha_b,
                                                            const StopToken& stop_token) {
  return shortenWithinHomotopy(
    polymap, alpha_a.projectedPath(), alpha_b.projectedPath(), stop_token);
}

HomotopyShorteningResult RaystarCore::shortenWithinHomotopy(const Polymap& polymap,
                                                            const std::vector<Point2d>& first,
                                                            const std::vector<Point2d>& second,
                                                            const StopToken& stop_token) {
  HomotopyShorteningResult invalid;
  if (stop_token.poll()) {
    invalid.status = HomotopyShorteningStatus::stopped;
    invalid.message = "UPS was canceled before composing the reference path";
    return invalid;
  }

  if (first.empty() || second.empty()) {
    invalid.status = HomotopyShorteningStatus::invalid_reference;
    invalid.message = "UPS requires two non-empty rooted references";
    return invalid;
  }
  if (first.front() != second.front()) {
    invalid.status = HomotopyShorteningStatus::invalid_reference;
    invalid.message = "UPS rooted references must have the same exact root point";
    return invalid;
  }

  size_t common_prefix = 0;
  while (common_prefix < first.size() && common_prefix < second.size() &&
         first[common_prefix] == second[common_prefix]) {
    if (stop_token.poll()) {
      invalid.status = HomotopyShorteningStatus::stopped;
      invalid.message = "UPS was canceled while removing the common reference prefix";
      return invalid;
    }
    ++common_prefix;
  }
  // The shared root check above guarantees at least one common point.
  const size_t branch_point = common_prefix - 1;
  std::vector<Point2d> reference;
  reference.reserve((first.size() - branch_point) + (second.size() - common_prefix));
  for (size_t index = first.size(); index > branch_point; --index)
    reference.push_back(first[index - 1]);
  for (size_t index = common_prefix; index < second.size(); ++index) {
    if (reference.empty() || reference.back() != second[index])
      reference.push_back(second[index]);
  }

  return polymap.shortenPathWithinHomotopy(reference, stop_token);
}

TransitionBatchResult RaystarCore::shortenReferenceTransitions(
  const Polymap& polymap,
  const std::vector<PathSolution>& configurations,
  const std::vector<ReferenceTransitionPair>& pairs,
  const StopToken& stop_token) {
  TransitionBatchResult result;
  result.transitions.reserve(pairs.size());
  if (stop_token.poll()) {
    result.status = TransitionBatchStatus::stopped;
    result.message = "UPS transition batch was canceled before validation";
    return result;
  }

  for (size_t index = 0; index < pairs.size(); ++index) {
    const auto& pair = pairs[index];
    if (pair.from_reference >= configurations.size() ||
        pair.to_reference >= configurations.size()) {
      result.status = TransitionBatchStatus::invalid_request;
      result.message = "UPS transition pair " + std::to_string(index) +
                       " contains an out-of-range configuration index";
      result.transitions.clear();
      return result;
    }
  }

  for (const auto& pair : pairs) {
    if (stop_token.poll()) {
      result.status = TransitionBatchStatus::stopped;
      result.message = "UPS transition batch was canceled after " +
                       std::to_string(result.transitions.size()) + " transitions";
      return result;
    }
    ReferenceTransitionResult transition;
    transition.pair = pair;
    transition.shortening = shortenWithinHomotopy(polymap,
                                                  configurations[pair.from_reference],
                                                  configurations[pair.to_reference],
                                                  stop_token);
    result.transitions.emplace_back(std::move(transition));
    if (result.transitions.back().shortening.status == HomotopyShorteningStatus::stopped) {
      result.status = TransitionBatchStatus::stopped;
      result.message = "UPS transition batch was canceled while shortening transition " +
                       std::to_string(result.transitions.size() - 1);
      return result;
    }
  }

  result.status = TransitionBatchStatus::success;
  result.message = "All requested UPS transitions were evaluated";
  return result;
}

PlanResult RaystarCore::plan(const GridMap& grid_map,
                             int start_x,
                             int start_y,
                             int goal_x,
                             int goal_y,
                             int K,
                             bool allow_self_crossing,
                             const PlanningLimits& limits) {
  return plan(grid_map,
              start_x,
              start_y,
              goal_x,
              goal_y,
              SearchObjective::topK(K),
              allow_self_crossing,
              limits);
}

PlanResult RaystarCore::plan(const GridMap& grid_map,
                             int start_x,
                             int start_y,
                             int goal_x,
                             int goal_y,
                             const SearchObjective& objective,
                             bool allow_self_crossing,
                             const PlanningLimits& limits) {
  std::string legacy_error;
  if (!validatePlanRequest(
        grid_map, start_x, start_y, goal_x, goal_y, objective, limits, legacy_error)) {
    std::vector<Node>().swap(N_);
    PlanResult result;
    result.outcome = PlanningOutcome::invalid_request;
    result.message = std::move(legacy_error);
    return result;
  }
  return plan(
    grid_map,
    PlanEndpoint(start_x, start_y, static_cast<double>(start_x), static_cast<double>(start_y)),
    PlanEndpoint(goal_x, goal_y, static_cast<double>(goal_x), static_cast<double>(goal_y)),
    objective,
    allow_self_crossing,
    limits);
}

PlanResult RaystarCore::plan(const GridMap& grid_map,
                             const Point2d& start,
                             const Point2d& goal,
                             int K,
                             bool allow_self_crossing,
                             const PlanningLimits& limits) {
  return plan(grid_map, start, goal, SearchObjective::topK(K), allow_self_crossing, limits);
}

PlanResult RaystarCore::plan(const GridMap& grid_map,
                             const Point2d& start,
                             const Point2d& goal,
                             const SearchObjective& objective,
                             bool allow_self_crossing,
                             const PlanningLimits& limits) {
  const auto endpoint_from_point = [&grid_map](const Point2d& point) {
    std::pair<int, int> cell = {-1, -1};
    if (std::isfinite(point.first) && std::isfinite(point.second)) {
      const double x = std::floor(point.first);
      const double y = std::floor(point.second);
      if (x >= static_cast<double>(std::numeric_limits<int>::lowest()) &&
          x <= static_cast<double>(std::numeric_limits<int>::max()) &&
          y >= static_cast<double>(std::numeric_limits<int>::lowest()) &&
          y <= static_cast<double>(std::numeric_limits<int>::max())) {
        cell = {static_cast<int>(x), static_cast<int>(y)};
      }
    }
    (void)grid_map;
    return PlanEndpoint(cell, point);
  };
  return plan(grid_map,
              endpoint_from_point(start),
              endpoint_from_point(goal),
              objective,
              allow_self_crossing,
              limits);
}

PlanResult RaystarCore::plan(const GridMap& grid_map,
                             const PlanEndpoint& start,
                             const PlanEndpoint& goal,
                             int K,
                             bool allow_self_crossing,
                             const PlanningLimits& limits) {
  return plan(grid_map, start, goal, SearchObjective::topK(K), allow_self_crossing, limits);
}

MultiGoalPlanResult RaystarCore::planToGoalsWithinCosts(const GridMap& grid_map,
                                                        const PlanEndpoint& start,
                                                        const std::vector<CostBoundedGoal>& goals,
                                                        bool allow_self_crossing,
                                                        const PlanningLimits& limits) try {
  resetSearchState();

  MultiGoalPlanResult result;
  result.goal_results.reserve(goals.size());
  for (const auto& goal : goals) {
    CostBoundedGoalResult goal_result;
    goal_result.endpoint = goal.endpoint;
    goal_result.requested_max_path_cost = goal.max_path_cost;
    result.goal_results.emplace_back(std::move(goal_result));
  }

  std::string validation_error;
  if (!validateMultiGoalPlanRequest(grid_map, start, goals, limits, validation_error)) {
    result.outcome = PlanningOutcome::invalid_request;
    result.message = validation_error;
    for (auto& goal_result : result.goal_results) {
      goal_result.outcome = PlanningOutcome::invalid_request;
      goal_result.message = validation_error;
    }
    return result;
  }

  for (auto& goal_result : result.goal_results) {
    const size_t reserve = std::min(limits.max_cost_bounded_paths, limits.max_path_points / 2);
    goal_result.path_solutions.reserve(reserve);
  }

  using PlanningClock = std::chrono::steady_clock;
  const auto request_start_time = PlanningClock::now();
  const bool timeout_enabled = limits.planning_timeout != std::chrono::milliseconds::max();
  PlanningLimitReached requested_stop_reason = PlanningLimitReached::none;
  const StopToken stop_token([&]() {
    if (limits.cancel_requested && limits.cancel_requested()) {
      requested_stop_reason = PlanningLimitReached::cancelled;
      return true;
    }
    if (timeout_enabled &&
        std::chrono::duration_cast<std::chrono::milliseconds>(
          PlanningClock::now() - request_start_time) >= limits.planning_timeout) {
      requested_stop_reason = PlanningLimitReached::timeout;
      return true;
    }
    return false;
  });

  std::vector<bool> finished(goals.size(), false);
  const auto sort_goal_paths = [&]() {
    for (auto& goal_result : result.goal_results) {
      std::stable_sort(goal_result.path_solutions.begin(),
                       goal_result.path_solutions.end(),
                       [](const PathSolution& lhs, const PathSolution& rhs) {
                         return lhs.path_cost_ < rhs.path_cost_;
                       });
    }
  };
  const auto total_path_count = [&]() {
    size_t count = 0;
    for (const auto& goal_result : result.goal_results) count += goal_result.path_solutions.size();
    return count;
  };
  const auto stop_for_limit = [&](PlanningLimitReached limit,
                                  std::shared_ptr<const Polymap> polymap,
                                  const std::optional<PlanningClock::time_point>& planner_start) {
    result.outcome = PlanningOutcome::limit_reached;
    result.limit_reached = limit;
    result.expanded_nodes = N_.size();
    result.message = planningLimitMessage(limit, limits, total_path_count());
    result.polymap = std::move(polymap);
    if (planner_start) {
      result.plan_time_ms =
        std::chrono::duration_cast<std::chrono::microseconds>(PlanningClock::now() - *planner_start)
          .count() /
        1000.0;
    }
    for (size_t i = 0; i < result.goal_results.size(); ++i) {
      if (finished[i])
        continue;
      auto& goal_result = result.goal_results[i];
      goal_result.success = !goal_result.path_solutions.empty();
      goal_result.outcome = PlanningOutcome::limit_reached;
      goal_result.limit_reached = limit;
      goal_result.message = planningLimitMessage(limit, limits, goal_result.path_solutions.size());
    }
    sort_goal_paths();
    return result;
  };
  const auto stop_for_request = [&](std::shared_ptr<const Polymap> polymap,
                                    const std::optional<PlanningClock::time_point>& planner_start) {
    return stop_for_limit(requested_stop_reason, std::move(polymap), planner_start);
  };

  if (stop_token.poll())
    return stop_for_request(nullptr, std::nullopt);

  GridMap work_map = grid_map;
  if (stop_token.poll())
    return stop_for_request(nullptr, std::nullopt);
  if (outlineMap(work_map.data, work_map.width, work_map.height, stop_token) ==
      OperationStatus::stopped) {
    return stop_for_request(nullptr, std::nullopt);
  }

  const int start_x = start.cell_.first;
  const int start_y = start.cell_.second;
  const double start_gx = start.position_.first;
  const double start_gy = start.position_.second;
  if (work_map.at(static_cast<unsigned int>(start_x), static_cast<unsigned int>(start_y)) != 0) {
    result.outcome = PlanningOutcome::invalid_request;
    result.message = "Invalid start: corresponding cell is occupied or on the map boundary";
    for (auto& goal_result : result.goal_results) {
      goal_result.outcome = PlanningOutcome::invalid_request;
      goal_result.message = result.message;
    }
    return result;
  }
  for (size_t i = 0; i < goals.size(); ++i) {
    const auto& cell = goals[i].endpoint.cell_;
    if (work_map.at(static_cast<unsigned int>(cell.first),
                    static_cast<unsigned int>(cell.second)) != 0) {
      result.outcome = PlanningOutcome::invalid_request;
      result.message = "Invalid goal[" + std::to_string(i) +
                       "]: corresponding cell is occupied or on the map boundary";
      for (auto& goal_result : result.goal_results) {
        goal_result.outcome = PlanningOutcome::invalid_request;
        goal_result.message = result.message;
      }
      return result;
    }
  }

  const auto map_start_time = PlanningClock::now();
  // Classify reachability once on the outlined occupancy map.  Unreachable
  // goals receive an individual no-path certificate; reachable goals alone are
  // passed to Polymap so one disconnected endpoint cannot discard the shared
  // geometry needed by the rest of the request.
  const size_t width = static_cast<size_t>(work_map.width);
  const size_t cell_count = width * static_cast<size_t>(work_map.height);
  std::vector<uint8_t> reachable(cell_count, 0);
  std::vector<int> flood_stack;
  flood_stack.reserve(std::min<size_t>(cell_count, 4096));
  const int start_index = start_x + start_y * work_map.width;
  reachable[static_cast<size_t>(start_index)] = 1;
  flood_stack.emplace_back(start_index);
  while (!flood_stack.empty()) {
    if (stop_token.poll())
      return stop_for_request(nullptr, std::nullopt);
    const int current = flood_stack.back();
    flood_stack.pop_back();
    const int x = current % work_map.width;
    const int y = current / work_map.width;
    const auto enqueue_free = [&](int neighbor) {
      const size_t index = static_cast<size_t>(neighbor);
      if (reachable[index] == 0 && work_map.data[index] == 0) {
        reachable[index] = 1;
        flood_stack.emplace_back(neighbor);
      }
    };
    if (x > 0)
      enqueue_free(current - 1);
    if (x + 1 < static_cast<int>(work_map.width))
      enqueue_free(current + 1);
    if (y > 0)
      enqueue_free(current - work_map.width);
    if (y + 1 < static_cast<int>(work_map.height))
      enqueue_free(current + work_map.width);
  }

  std::vector<size_t> reachable_goal_indices;
  std::vector<PolymapEndpoint> polymap_goals;
  reachable_goal_indices.reserve(goals.size());
  polymap_goals.reserve(goals.size());
  for (size_t i = 0; i < goals.size(); ++i) {
    const auto& endpoint = goals[i].endpoint;
    const size_t index = static_cast<size_t>(endpoint.cell_.second) * width +
                         static_cast<size_t>(endpoint.cell_.first);
    if (reachable[index] == 0) {
      auto& goal_result = result.goal_results[i];
      goal_result.outcome = PlanningOutcome::no_path;
      goal_result.completion = PlanningCompletion::frontier_exhausted;
      goal_result.message = "No path exists between start and this goal";
      finished[i] = true;
      continue;
    }
    reachable_goal_indices.emplace_back(i);
    polymap_goals.push_back({endpoint.cell_.first, endpoint.cell_.second, endpoint.position_});
  }
  // Do not overlap the Core classification scratch buffers with Polymap's own
  // flood-fill mask. The documented map working-set estimate budgets one such
  // reachability phase at a time.
  std::vector<uint8_t>().swap(reachable);
  std::vector<int>().swap(flood_stack);

  if (reachable_goal_indices.empty()) {
    result.map_time_ms =
      std::chrono::duration_cast<std::chrono::microseconds>(PlanningClock::now() - map_start_time)
        .count() /
      1000.0;
    result.outcome = PlanningOutcome::no_path;
    result.message = "No requested goal is reachable from start";
    return result;
  }

  auto map_build =
    Polymap::create(work_map, start_x, start_y, start.position_, polymap_goals, stop_token, limits);
  const auto map_end_time = PlanningClock::now();
  result.map_time_ms =
    std::chrono::duration_cast<std::chrono::microseconds>(map_end_time - map_start_time).count() /
    1000.0;
  if (map_build.status == PolymapCreateStatus::stopped || stop_token.poll())
    return stop_for_request(nullptr, std::nullopt);
  if (map_build.status == PolymapCreateStatus::no_path) {
    result.outcome = PlanningOutcome::failed;
    result.message = "Polymap reachability disagrees with the shared reachability classification";
    for (const size_t i : reachable_goal_indices) {
      result.goal_results[i].outcome = PlanningOutcome::failed;
      result.goal_results[i].message = result.message;
    }
    return result;
  }
  if (map_build.status == PolymapCreateStatus::failure || !map_build) {
    constexpr const char* invalid_start_prefix = "Invalid start position: ";
    constexpr const char* invalid_goal_prefix = "Invalid goal position: ";
    if (map_build.error.rfind(invalid_start_prefix, 0) == 0 ||
        map_build.error.rfind(invalid_goal_prefix, 0) == 0) {
      result.outcome = PlanningOutcome::invalid_request;
      result.message = "Invalid multi-goal endpoint: " + map_build.error;
    } else {
      result.outcome = PlanningOutcome::failed;
      result.message = "Map geometry construction failed";
      if (!map_build.error.empty())
        result.message += ": " + map_build.error;
    }
    for (const size_t i : reachable_goal_indices) {
      result.goal_results[i].outcome = result.outcome;
      result.goal_results[i].message = result.message;
    }
    return result;
  }

  auto theMap = std::make_shared<Polymap>(std::move(*map_build.value));
  const auto planner_start_time = PlanningClock::now();
  const std::optional<PlanningClock::time_point> planner_start = planner_start_time;

  const auto fail_planning = [&](const std::string& message,
                                 PlanningOutcome outcome = PlanningOutcome::failed) {
    result.outcome = outcome;
    result.expanded_nodes = N_.size();
    result.message = message;
    result.plan_time_ms = std::chrono::duration_cast<std::chrono::microseconds>(
                            PlanningClock::now() - planner_start_time)
                            .count() /
                          1000.0;
    result.polymap = theMap;
    for (size_t i = 0; i < result.goal_results.size(); ++i) {
      if (finished[i])
        continue;
      auto& goal_result = result.goal_results[i];
      goal_result.success = false;
      goal_result.path_solutions.clear();
      goal_result.outcome = outcome;
      goal_result.message = message;
    }
    resetSearchState();
    return result;
  };

  std::vector<bool> active(goals.size(), false);
  for (const size_t i : reachable_goal_indices) active[i] = true;
  size_t active_count = reachable_goal_indices.size();
  size_t accumulated_path_points = 0;

  const auto goal_lower_bound =
    [&](const Candidate& candidate, size_t goal_index, double& lower_bound) {
      ConservativeBinary64PathLength certificate(maximum_length_refinement_precision_);
      return buildCandidateLengthCertificate(start.position_,
                                             N_,
                                             candidate,
                                             goals[goal_index].endpoint.position_,
                                             certificate,
                                             &stop_token) &&
             certificate.lowerBound(lower_bound);
    };
  const auto signed_bound_slack = [](double lower_bound, double inclusive_bound) {
    const double rounded = lower_bound - inclusive_bound;
    // The subtraction is used only as a heap key. Clamp its sign from the
    // exact binary64 operand comparison, so no rounding mode can turn an
    // eligible candidate into a positive (prunable) key or vice versa.
    if (lower_bound <= inclusive_bound)
      return std::min(rounded, 0.0);
    return rounded > 0.0 ? rounded : std::numeric_limits<double>::denorm_min();
  };
  const auto shared_candidate_key = [&](const Candidate& candidate, double& key) {
    key = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < goals.size(); ++i) {
      if (!active[i])
        continue;
      double lower_bound = std::numeric_limits<double>::quiet_NaN();
      if (!goal_lower_bound(candidate, i, lower_bound))
        return false;
      key = std::min(key, signed_bound_slack(lower_bound, goals[i].max_path_cost));
    }
    return true;
  };
  const auto minimum_goal_distance = [&](const Point2d& source) {
    double distance = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < goals.size(); ++i) {
      if (!active[i])
        continue;
      const auto& goal_position = goals[i].endpoint.position_;
      distance = std::min(
        distance,
        std::hypot(source.first - goal_position.first, source.second - goal_position.second));
    }
    return distance;
  };
  const auto comp = [](const Candidate& a, const Candidate& b) { return a.Fcost_ > b.Fcost_; };
  const auto complete_goal = [&](size_t i, PlanningCompletion completion) {
    auto& goal_result = result.goal_results[i];
    goal_result.success = !goal_result.path_solutions.empty();
    goal_result.outcome =
      goal_result.success ? PlanningOutcome::complete : PlanningOutcome::no_path;
    goal_result.completion = completion;
    goal_result.message = goal_result.success
                            ? "Cost-bounded enumeration complete; found " +
                                std::to_string(goal_result.path_solutions.size()) + " path(s)"
                            : "Cost-bounded enumeration complete; no path found";
    active[i] = false;
    finished[i] = true;
    --active_count;
  };

  std::vector<Candidate> queue;
  Candidate root_candidate(-1, -1, 0.0);
  if (!shared_candidate_key(root_candidate, root_candidate.Fcost_)) {
    if (stop_token.poll())
      return stop_for_request(theMap, planner_start);
    return fail_planning("Could not certify the root frontier lower bound");
  }
  queue.emplace_back(root_candidate);

  const auto rebuild_queue = [&]() {
    for (auto& candidate : queue) {
      if (!shared_candidate_key(candidate, candidate.Fcost_))
        return false;
    }
    std::make_heap(queue.begin(), queue.end(), comp);
    return true;
  };

  while (!queue.empty() && active_count > 0) {
    if (stop_token.poll())
      return stop_for_request(theMap, planner_start);
    if (queue.front().Fcost_ > 0.0) {
      for (size_t i = 0; i < goals.size(); ++i) {
        if (active[i])
          complete_goal(i, PlanningCompletion::cost_bound_exhausted);
      }
      break;
    }
    if (N_.size() >= limits.max_nodes)
      return stop_for_limit(PlanningLimitReached::max_nodes, theMap, planner_start);

    Candidate best_candidate = queue.front();
    std::pop_heap(queue.begin(), queue.end(), comp);
    queue.pop_back();
    if (stop_token.poll())
      return stop_for_request(theMap, planner_start);

    const int parent_index = best_candidate.Nindex_;
    const int child_index = best_candidate.Cindex_;
    std::optional<Node> pending_node;
    if (parent_index == -1) {
      VisibilityRegion visibility;
      std::string visibility_error;
      const Point2d start_position = {start_gx, start_gy};
      const auto visibility_status =
        theMap->getRootVisibilityRegion(start_position, visibility, stop_token, &visibility_error);
      if (visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request(theMap, planner_start);
      if (visibility_status == OperationStatus::failure)
        return fail_planning("Root visibility calculation failed: " + visibility_error);
      if (N_.size() >= limits.max_nodes)
        return stop_for_limit(PlanningLimitReached::max_nodes, theMap, planner_start);

      Node new_node(theMap.get(),
                    0,
                    static_cast<double>(start_x),
                    static_cast<double>(start_y),
                    0.0,
                    minimum_goal_distance(start_position),
                    visibility,
                    stop_token);
      new_node.is_continuous_root_ = true;
      new_node.local_shortest_path_.clear();
      const auto full_visibility_status = new_node.setFullVisibilityRegion(visibility, stop_token);
      if (full_visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request(theMap, planner_start);
      if (full_visibility_status == OperationStatus::failure)
        return fail_planning("Root full visibility projection failed");

      std::string child_error;
      const auto child_status =
        generateChildrenFromSource(theMap.get(),
                                   new_node.Nindex_,
                                   start_position,
                                   exact_geometry::Point(start_gx, start_gy),
                                   new_node.start_angle_,
                                   new_node.end_angle_,
                                   new_node.Gcost_,
                                   new_node.visibility_region_,
                                   new_node.C_,
                                   stop_token,
                                   &child_error);
      if (child_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request(theMap, planner_start);
      if (child_status == OperationStatus::failure)
        return fail_planning("Root child generation failed: " + child_error);
      pending_node.emplace(std::move(new_node));
    } else {
      if (parent_index < 0 || parent_index >= static_cast<int>(N_.size()))
        return fail_planning("Planner candidate has an invalid parent index");
      if (child_index < 0 ||
          child_index >= static_cast<int>(N_[static_cast<size_t>(parent_index)].C_.size())) {
        return fail_planning("Planner candidate has an invalid child index");
      }

      const auto new_source_point = N_[parent_index].C_[child_index].c_;
      const int new_node_index = static_cast<int>(N_.size());
      bool path_self_crosses = false;
      if (!allow_self_crossing) {
        const auto crossing_status =
          hybridPathSelfCrosses(start.position_,
                                N_[parent_index].local_shortest_path_,
                                {static_cast<double>(new_source_point.first),
                                 static_cast<double>(new_source_point.second)},
                                path_self_crosses,
                                stop_token);
        if (crossing_status == OperationStatus::stopped || stop_token.poll())
          return stop_for_request(theMap, planner_start);
        if (crossing_status == OperationStatus::failure)
          return fail_planning("Path self-crossing validation failed");
      }
      if (path_self_crosses)
        continue;

      VisibilityRegion visibility;
      VisibilityRegion full_visibility;
      std::string visibility_error;
      const auto visibility_status = getScopedVisibilityRegion(
        *theMap, best_candidate, visibility, full_visibility, stop_token, visibility_error);
      if (visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request(theMap, planner_start);
      if (visibility_status == OperationStatus::failure)
        return fail_planning("Scoped visibility calculation failed: " + visibility_error);

      Node new_node(theMap.get(),
                    new_node_index,
                    new_source_point.first,
                    new_source_point.second,
                    N_[parent_index].C_[child_index].c_gcost_,
                    N_[parent_index].C_[child_index].c_hcost_,
                    parent_index,
                    visibility,
                    stop_token);
      const auto full_visibility_status =
        new_node.setFullVisibilityRegion(full_visibility, stop_token);
      if (full_visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request(theMap, planner_start);
      if (full_visibility_status == OperationStatus::failure)
        return fail_planning("Scoped full visibility projection failed");

      new_node.local_shortest_path_.reserve(N_[parent_index].local_shortest_path_.size() + 1);
      new_node.local_shortest_path_ = N_[parent_index].local_shortest_path_;
      new_node.local_shortest_path_.emplace_back(new_source_point);
      new_node.path_node_index_ = N_[parent_index].path_node_index_;
      new_node.path_node_index_.emplace_back(new_node_index);
      new_node.start_angle_ = N_[parent_index].C_[child_index].start_angle_;
      new_node.end_angle_ = N_[parent_index].C_[child_index].end_angle_;

      std::string child_error;
      const auto child_status = new_node.generateChild(theMap.get(), stop_token, &child_error);
      if (child_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request(theMap, planner_start);
      if (child_status == OperationStatus::failure)
        return fail_planning("Scoped child generation failed: " + child_error);
      pending_node.emplace(std::move(new_node));
    }

    if (!pending_node)
      return fail_planning("Planner failed to construct the selected node");
    for (auto& child : pending_node->C_)
      child.c_hcost_ = minimum_goal_distance(child.c_endpoint_.position);

    N_.emplace_back(std::move(*pending_node));
    for (const auto& child : N_.back().C_) {
      Candidate candidate(child.Nindex_, child.Cindex_, 0.0);
      if (!shared_candidate_key(candidate, candidate.Fcost_)) {
        if (stop_token.poll())
          return stop_for_request(theMap, planner_start);
        return fail_planning("Could not certify a frontier path lower bound");
      }
      queue.emplace_back(candidate);
      std::push_heap(queue.begin(), queue.end(), comp);
    }

    const auto point_location_mode = visibilityBoundaryModeForNode(*theMap, N_.back());
    const exact_geometry::Point exact_node_source =
      N_.back().is_continuous_root_
        ? exact_geometry::Point(start_gx, start_gy)
        : exact_geometry::Point(N_.back().seed_.first, N_.back().seed_.second);

    bool active_set_changed = false;
    for (size_t i = 0; i < goals.size(); ++i) {
      if (!active[i])
        continue;
      if (stop_token.poll())
        return stop_for_request(theMap, planner_start);
      const auto& goal_position = goals[i].endpoint.position_;
      std::pair<bool, bool> location = {false, false};
      const auto point_location_status = detail::classifyPointInVisibilityBoundary(
        N_.back().visibility_region_,
        exact_geometry::Point(goal_position.first, goal_position.second),
        exact_node_source,
        point_location_mode,
        location,
        stop_token);
      if (point_location_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request(theMap, planner_start);
      if (point_location_status == OperationStatus::failure)
        return fail_planning("Goal[" + std::to_string(i) + "] point-location validation failed");
      if (!location.first && !location.second)
        continue;

      bool goal_segment_crosses = false;
      if (!allow_self_crossing) {
        const auto crossing_status = hybridPathSelfCrosses(start.position_,
                                                           N_.back().local_shortest_path_,
                                                           goal_position,
                                                           goal_segment_crosses,
                                                           stop_token);
        if (crossing_status == OperationStatus::stopped || stop_token.poll())
          return stop_for_request(theMap, planner_start);
        if (crossing_status == OperationStatus::failure)
          return fail_planning("Goal path self-crossing validation failed");
      }
      if (goal_segment_crosses)
        continue;

      ConservativeBinary64PathLength path_certificate(maximum_length_refinement_precision_);
      if (!buildPolylineLengthCertificate(start.position_,
                                          N_.back().local_shortest_path_,
                                          goal_position,
                                          path_certificate,
                                          &stop_token)) {
        if (stop_token.poll())
          return stop_for_request(theMap, planner_start);
        return fail_planning("Could not construct a goal path length certificate");
      }
      ConservativeBinary64PathLength::Comparison bound_comparison =
        ConservativeBinary64PathLength::Comparison::equal;
      const auto certificate_stop = [&]() { return stop_token.poll(); };
      if (!path_certificate.compareTo(goals[i].max_path_cost, bound_comparison, certificate_stop)) {
        if (stop_token.poll())
          return stop_for_request(theMap, planner_start);
        return fail_planning("Could not resolve a goal path against its inclusive cost bound");
      }
      if (bound_comparison == ConservativeBinary64PathLength::Comparison::greater)
        continue;
      double path_cost = std::numeric_limits<double>::quiet_NaN();
      if (!path_certificate.upperBound(path_cost, certificate_stop) ||
          path_cost > goals[i].max_path_cost) {
        if (stop_token.poll())
          return stop_for_request(theMap, planner_start);
        return fail_planning("Goal path upper certificate disagrees with exact bound comparison");
      }

      if (goals[i].path_admission) {
        const BoundedPathView path_view{
          start.position_, N_.back().local_shortest_path_, goal_position, path_cost};
        const auto admission = goals[i].path_admission(path_view, stop_token);
        // A callback may observe and latch the cooperative stop while still
        // returning any enum value. Stop always outranks subsequent path or
        // point resource admission.
        if (stop_token.poll())
          return stop_for_request(theMap, planner_start);
        if (admission == BoundedPathAdmission::reject)
          continue;
        if (admission != BoundedPathAdmission::accept) {
          return fail_planning("Goal[" + std::to_string(i) + "] bounded path admission failed");
        }
      }

      // The cap is charged only to candidates accepted by the stronger
      // adapter certificate. Continue expanding after retaining exactly the
      // cap so rejected search-superset shells can still be exhausted; stop
      // only when another eligible path would exceed it.
      auto& goal_result = result.goal_results[i];
      if (goal_result.path_solutions.size() >= limits.max_cost_bounded_paths) {
        goal_result.success = true;
        goal_result.outcome = PlanningOutcome::limit_reached;
        goal_result.limit_reached = PlanningLimitReached::max_paths;
        goal_result.message = planningLimitMessage(
          PlanningLimitReached::max_paths, limits, goal_result.path_solutions.size());
        active[i] = false;
        finished[i] = true;
        --active_count;
        active_set_changed = true;
        continue;
      }

      const size_t path_points = N_.back().local_shortest_path_.size() + 2;
      if (path_points > limits.max_path_points ||
          accumulated_path_points > limits.max_path_points - path_points) {
        return stop_for_limit(PlanningLimitReached::max_path_points, theMap, planner_start);
      }
      goal_result.path_solutions.emplace_back(start.position_,
                                              N_.back().local_shortest_path_,
                                              goal_position,
                                              path_cost,
                                              N_.back().path_node_index_);
      accumulated_path_points += path_points;
    }

    // Retaining exactly the cap is not itself a limit: the remaining padded
    // frontier may contain only adapter-rejected shell candidates. Preserve
    // the original per-goal proof shortcut when its own frontier is already
    // exhausted, otherwise keep searching and report max_paths only if a
    // subsequent accepted solution actually attempts to exceed the cap.
    for (size_t i = 0; i < goals.size(); ++i) {
      if (!active[i] ||
          result.goal_results[i].path_solutions.size() < limits.max_cost_bounded_paths) {
        continue;
      }
      double frontier_lower_bound = std::numeric_limits<double>::infinity();
      for (const auto& candidate : queue) {
        double candidate_lower_bound = std::numeric_limits<double>::quiet_NaN();
        if (!goal_lower_bound(candidate, i, candidate_lower_bound)) {
          if (stop_token.poll())
            return stop_for_request(theMap, planner_start);
          return fail_planning("Could not certify a remaining frontier path lower bound");
        }
        frontier_lower_bound = std::min(frontier_lower_bound, candidate_lower_bound);
      }
      if (queue.empty()) {
        complete_goal(i, PlanningCompletion::frontier_exhausted);
        active_set_changed = true;
      } else if (frontier_lower_bound > goals[i].max_path_cost) {
        complete_goal(i, PlanningCompletion::cost_bound_exhausted);
        active_set_changed = true;
      }
    }

    if (active_set_changed && active_count > 0 && !rebuild_queue()) {
      if (stop_token.poll())
        return stop_for_request(theMap, planner_start);
      return fail_planning("Could not rebuild certified shared-frontier keys");
    }
  }

  if (stop_token.poll())
    return stop_for_request(theMap, planner_start);
  if (queue.empty()) {
    for (size_t i = 0; i < goals.size(); ++i) {
      if (active[i])
        complete_goal(i, PlanningCompletion::frontier_exhausted);
    }
  }

  sort_goal_paths();

  result.expanded_nodes = N_.size();
  result.plan_time_ms =
    std::chrono::duration_cast<std::chrono::microseconds>(PlanningClock::now() - planner_start_time)
      .count() /
    1000.0;
  result.polymap = theMap;
  bool any_path = false;
  bool any_limit = false;
  for (const auto& goal_result : result.goal_results) {
    any_path = any_path || !goal_result.path_solutions.empty();
    any_limit = any_limit || goal_result.outcome == PlanningOutcome::limit_reached;
  }
  if (any_limit) {
    result.outcome = PlanningOutcome::limit_reached;
    result.limit_reached = PlanningLimitReached::max_paths;
    result.message = "One or more goals reached max_cost_bounded_paths; other goals were completed";
  } else if (any_path) {
    result.outcome = PlanningOutcome::complete;
    result.message = "Multi-goal cost-bounded enumeration complete";
  } else {
    result.outcome = PlanningOutcome::no_path;
    result.message = "Multi-goal cost-bounded enumeration complete; no path found";
  }
  return result;
} catch (const std::exception& exception) {
  resetSearchState();
  MultiGoalPlanResult result;
  result.message = boundedPlanningExceptionMessage(exception.what());
  return result;
} catch (...) {
  resetSearchState();
  MultiGoalPlanResult result;
  result.message = "Planning failed with an unknown exception";
  return result;
}

PlanResult RaystarCore::plan(const GridMap& grid_map,
                             const PlanEndpoint& start,
                             const PlanEndpoint& goal,
                             const SearchObjective& objective,
                             bool allow_self_crossing,
                             const PlanningLimits& limits) try {
  if (objective.mode == SearchMode::all_within_cost) {
    std::string validation_error;
    if (!validatePlanRequest(grid_map, start, goal, objective, limits, validation_error)) {
      resetSearchState();
      PlanResult invalid_result;
      invalid_result.outcome = PlanningOutcome::invalid_request;
      invalid_result.message = std::move(validation_error);
      return invalid_result;
    }
    auto multi_result =
      planToGoalsWithinCosts(grid_map,
                             start,
                             {{goal, objective.max_path_cost, objective.path_admission}},
                             allow_self_crossing,
                             limits);
    PlanResult single_result;
    if (!multi_result.goal_results.empty()) {
      auto& goal_result = multi_result.goal_results.front();
      single_result.success = goal_result.success;
      single_result.message = std::move(goal_result.message);
      single_result.path_solutions = std::move(goal_result.path_solutions);
      single_result.outcome = goal_result.outcome;
      single_result.limit_reached = goal_result.limit_reached;
      single_result.completion = goal_result.completion;
    } else {
      single_result.message = multi_result.message;
      single_result.outcome = multi_result.outcome;
      single_result.limit_reached = multi_result.limit_reached;
    }
    single_result.expanded_nodes = multi_result.expanded_nodes;
    single_result.map_time_ms = multi_result.map_time_ms;
    single_result.plan_time_ms = multi_result.plan_time_ms;
    single_result.polymap = std::move(multi_result.polymap);
    return single_result;
  }

  // Per-plan state must never survive into a new request, including invalid
  // requests and maps for which Polymap reports that no solution exists.
  // Releasing capacity here prevents a previously large request from keeping
  // its entire debug/search allocation resident while the next map is copied.
  std::vector<Node>().swap(N_);

  PlanResult result;
  std::string validation_error;
  if (!validatePlanRequest(grid_map, start, goal, objective, limits, validation_error)) {
    result.outcome = PlanningOutcome::invalid_request;
    result.message = std::move(validation_error);
    return result;
  }

  // Every returned path consumes at least the two continuous endpoints.  A
  // bounded reserve avoids repeated geometric growth while still respecting
  // both K and the response-wide path-point budget.
  // Cost-bounded requests returned through the shared-tree implementation
  // above. From this point onward there is exactly one top-K search path; do
  // not retain a second, numerically divergent bounded implementation.
  const size_t requested_reserve = static_cast<size_t>(objective.k);
  const size_t path_reserve = std::min(requested_reserve, limits.max_path_points / 2);
  result.path_solutions.reserve(path_reserve);

  const int start_x = start.cell_.first;
  const int start_y = start.cell_.second;
  const int goal_x = goal.cell_.first;
  const int goal_y = goal.cell_.second;
  const double start_gx = start.position_.first;
  const double start_gy = start.position_.second;
  const double goal_gx = goal.position_.first;
  const double goal_gy = goal.position_.second;

  using PlanningClock = std::chrono::steady_clock;
  const auto request_start_time = PlanningClock::now();
  const bool timeout_enabled = limits.planning_timeout != std::chrono::milliseconds::max();
  PlanningLimitReached requested_stop_reason = PlanningLimitReached::none;
  const StopToken stop_token([&]() {
    if (limits.cancel_requested && limits.cancel_requested()) {
      requested_stop_reason = PlanningLimitReached::cancelled;
      return true;
    }
    if (timeout_enabled &&
        std::chrono::duration_cast<std::chrono::milliseconds>(
          PlanningClock::now() - request_start_time) >= limits.planning_timeout) {
      requested_stop_reason = PlanningLimitReached::timeout;
      return true;
    }
    return false;
  });

  const auto stop_before_planner = [&](std::shared_ptr<const Polymap> polymap) {
    result.outcome = PlanningOutcome::limit_reached;
    result.limit_reached = requested_stop_reason;
    result.message =
      planningLimitMessage(result.limit_reached, limits, result.path_solutions.size());
    result.polymap = std::move(polymap);
    return result;
  };

  if (stop_token.poll())
    return stop_before_planner(nullptr);

  GridMap work_map = grid_map;
  if (stop_token.poll())
    return stop_before_planner(nullptr);
  if (outlineMap(work_map.data, work_map.width, work_map.height, stop_token) ==
      OperationStatus::stopped) {
    return stop_before_planner(nullptr);
  }

  // Outlining reserves the map border as a geometric boundary.  Reject an
  // endpoint whose cell was consumed by that operation instead of silently
  // clearing it (the former behaviour could turn an occupied start into free
  // space and made a boundary point appear valid).
  if (work_map.at(static_cast<unsigned int>(start_x), static_cast<unsigned int>(start_y)) != 0) {
    result.outcome = PlanningOutcome::invalid_request;
    result.message = "Invalid start: corresponding cell is occupied or on the map boundary";
    return result;
  }
  if (work_map.at(static_cast<unsigned int>(goal_x), static_cast<unsigned int>(goal_y)) != 0) {
    result.outcome = PlanningOutcome::invalid_request;
    result.message = "Invalid goal: corresponding cell is occupied or on the map boundary";
    return result;
  }

  const auto map_start_time = PlanningClock::now();
  auto map_build = Polymap::create(work_map,
                                   start_x,
                                   start_y,
                                   goal_x,
                                   goal_y,
                                   start.position_,
                                   goal.position_,
                                   stop_token,
                                   limits);
  const auto map_end_time = PlanningClock::now();
  result.map_time_ms =
    std::chrono::duration_cast<std::chrono::microseconds>(map_end_time - map_start_time).count() /
    1000.0;

  // Poll after construction as well: a single CGAL primitive cannot be
  // pre-empted, so the deadline may have elapsed inside that call.
  if (map_build.status == PolymapCreateStatus::stopped || stop_token.poll())
    return stop_before_planner(nullptr);

  if (map_build.status == PolymapCreateStatus::no_path) {
    result.success = false;
    result.outcome = PlanningOutcome::no_path;
    result.completion = PlanningCompletion::frontier_exhausted;
    result.message = "No path exists between start and goal";
    return result;
  }

  if (map_build.status == PolymapCreateStatus::failure || !map_build) {
    result.success = false;
    // Endpoint validation in the continuous Polymap constructor is an input
    // contract failure, not a generic CDT failure.  Preserve the precise
    // start/goal diagnostic so callers can correct the request directly.
    constexpr const char* invalid_start_prefix = "Invalid start position: ";
    constexpr const char* invalid_goal_prefix = "Invalid goal position: ";
    if (map_build.error.rfind(invalid_start_prefix, 0) == 0) {
      result.outcome = PlanningOutcome::invalid_request;
      result.message = "Invalid start: point must be a strict free-space interior; " +
                       map_build.error.substr(std::string(invalid_start_prefix).size());
    } else if (map_build.error.rfind(invalid_goal_prefix, 0) == 0) {
      result.outcome = PlanningOutcome::invalid_request;
      result.message = "Invalid goal: point must be a strict free-space interior; " +
                       map_build.error.substr(std::string(invalid_goal_prefix).size());
    } else {
      result.outcome = PlanningOutcome::failed;
      result.message = "Map geometry construction failed";
    }
    if (!map_build.error.empty()) {
      if (result.message == "Map geometry construction failed")
        result.message += ": " + map_build.error;
    }
    return result;
  }

  auto theMap = std::make_shared<Polymap>(std::move(*map_build.value));

  const auto planner_start_time = PlanningClock::now();

  const auto fail_planning = [&](const std::string& message,
                                 PlanningOutcome outcome = PlanningOutcome::failed) {
    const auto failure_time = PlanningClock::now();
    result.success = false;
    result.outcome = outcome;
    result.expanded_nodes = N_.size();
    result.message = message;
    result.path_solutions.clear();
    result.plan_time_ms =
      std::chrono::duration_cast<std::chrono::microseconds>(failure_time - planner_start_time)
        .count() /
      1000.0;
    result.polymap = theMap;
    resetSearchState();
    return result;
  };

  std::string endpoint_error;
  const auto start_interior_status =
    theMap->validateFreeSpaceInterior(start.position_, stop_token, &endpoint_error);
  if (start_interior_status == OperationStatus::stopped || stop_token.poll())
    return stop_before_planner(theMap);
  if (start_interior_status == OperationStatus::failure) {
    return fail_planning(
      "Invalid start: point must be a strict free-space interior; " + endpoint_error,
      PlanningOutcome::invalid_request);
  }
  const auto goal_interior_status =
    theMap->validateFreeSpaceInterior(goal.position_, stop_token, &endpoint_error);
  if (goal_interior_status == OperationStatus::stopped || stop_token.poll())
    return stop_before_planner(theMap);
  if (goal_interior_status == OperationStatus::failure) {
    return fail_planning(
      "Invalid goal: point must be a strict free-space interior; " + endpoint_error,
      PlanningOutcome::invalid_request);
  }

  const auto stop_for_limit = [&](PlanningLimitReached limit) {
    const auto stop_time = PlanningClock::now();
    result.outcome = PlanningOutcome::limit_reached;
    result.limit_reached = limit;
    result.expanded_nodes = N_.size();
    result.success = !result.path_solutions.empty();
    result.message = planningLimitMessage(limit, limits, result.path_solutions.size());
    result.plan_time_ms =
      std::chrono::duration_cast<std::chrono::microseconds>(stop_time - planner_start_time)
        .count() /
      1000.0;
    result.polymap = theMap;
    return result;
  };

  const auto stop_for_request = [&]() { return stop_for_limit(requested_stop_reason); };

  size_t accumulated_path_points = 0;

  auto comp = [](const Candidate& a, const Candidate& b) { return a.Fcost_ > b.Fcost_; };

  std::vector<Candidate> queue;
  queue.emplace_back(Candidate(-1, -1, std::hypot(start_gx - goal_gx, start_gy - goal_gy)));

  while (!queue.empty()) {
    if (stop_token.poll())
      return stop_for_request();
    if (N_.size() >= limits.max_nodes)
      return stop_for_limit(PlanningLimitReached::max_nodes);

    Candidate best_candidate = queue[0];
    std::pop_heap(queue.begin(), queue.end(), comp);
    queue.pop_back();
    if (stop_token.poll())
      return stop_for_request();

    int parent_index = best_candidate.Nindex_;
    int child_index = best_candidate.Cindex_;
    std::optional<Node> pending_node;

    if (parent_index == -1) {
      VisibilityRegion Vtemp;
      std::string visibility_error;
      const Point2d start_position = {start_gx, start_gy};
      const auto visibility_status =
        theMap->getRootVisibilityRegion(start_position, Vtemp, stop_token, &visibility_error);
      if (visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (visibility_status == OperationStatus::failure) {
        return fail_planning("Root visibility calculation failed: " + visibility_error);
      }
      if (N_.size() >= limits.max_nodes)
        return stop_for_limit(PlanningLimitReached::max_nodes);

      Node new_node(theMap.get(),
                    0,
                    static_cast<double>(start_x),
                    static_cast<double>(start_y),
                    0.0,
                    best_candidate.Fcost_,
                    Vtemp,
                    stop_token);
      new_node.is_continuous_root_ = true;
      // The root holder remains in N_[0] for the legacy debug/index contract,
      // but its integer seed is not a geometric waypoint.
      new_node.local_shortest_path_.clear();
      if (stop_token.poll())
        return stop_for_request();
      const auto full_visibility_status = new_node.setFullVisibilityRegion(Vtemp, stop_token);
      if (full_visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (full_visibility_status == OperationStatus::failure)
        return fail_planning("Root full visibility projection failed");

      std::string child_error;
      const auto child_status =
        generateChildrenFromSource(theMap.get(),
                                   new_node.Nindex_,
                                   start_position,
                                   exact_geometry::Point(start_gx, start_gy),
                                   new_node.start_angle_,
                                   new_node.end_angle_,
                                   new_node.Gcost_,
                                   new_node.visibility_region_,
                                   new_node.C_,
                                   stop_token,
                                   &child_error);
      if (child_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (child_status == OperationStatus::failure)
        return fail_planning("Root child generation failed: " + child_error);
      if (N_.size() >= limits.max_nodes)
        return stop_for_limit(PlanningLimitReached::max_nodes);
      pending_node.emplace(std::move(new_node));
    } else {
      if (parent_index < 0 || parent_index >= static_cast<int>(N_.size()))
        return fail_planning("Planner candidate has an invalid parent index");
      if (child_index < 0 ||
          child_index >= static_cast<int>(N_[static_cast<size_t>(parent_index)].C_.size()))
        return fail_planning("Planner candidate has an invalid child index");

      auto new_source_point = N_[parent_index].C_[child_index].c_;
      int new_node_index = static_cast<int>(N_.size());

      bool path_self_crosses = false;
      if (!allow_self_crossing) {
        const auto crossing_status =
          hybridPathSelfCrosses({start_gx, start_gy},
                                N_[parent_index].local_shortest_path_,
                                {static_cast<double>(new_source_point.first),
                                 static_cast<double>(new_source_point.second)},
                                path_self_crosses,
                                stop_token);
        if (crossing_status == OperationStatus::stopped || stop_token.poll())
          return stop_for_request();
        if (crossing_status == OperationStatus::failure)
          return fail_planning("Path self-crossing validation failed");
      }
      if (path_self_crosses)
        continue;

      VisibilityRegion Vtemp;
      VisibilityRegion fullVtemp;
      std::string visibility_error;
      const auto visibility_status = getScopedVisibilityRegion(
        *theMap, best_candidate, Vtemp, fullVtemp, stop_token, visibility_error);
      if (visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (visibility_status == OperationStatus::failure) {
        return fail_planning("Scoped visibility calculation failed: " + visibility_error);
      }
      if (N_.size() >= limits.max_nodes)
        return stop_for_limit(PlanningLimitReached::max_nodes);

      Node new_node(theMap.get(),
                    new_node_index,
                    new_source_point.first,
                    new_source_point.second,
                    N_[parent_index].C_[child_index].c_gcost_,
                    N_[parent_index].C_[child_index].c_hcost_,
                    parent_index,
                    Vtemp,
                    stop_token);
      if (stop_token.poll())
        return stop_for_request();
      const auto full_visibility_status = new_node.setFullVisibilityRegion(fullVtemp, stop_token);
      if (full_visibility_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (full_visibility_status == OperationStatus::failure)
        return fail_planning("Scoped full visibility projection failed");

      new_node.local_shortest_path_.clear();
      new_node.local_shortest_path_.reserve(N_[parent_index].local_shortest_path_.size() + 1);
      for (const auto& waypoint : N_[parent_index].local_shortest_path_) {
        if (stop_token.poll())
          return stop_for_request();
        new_node.local_shortest_path_.emplace_back(waypoint);
      }
      new_node.local_shortest_path_.emplace_back(new_source_point);
      new_node.path_node_index_ = N_[parent_index].path_node_index_;
      new_node.path_node_index_.emplace_back(new_node_index);
      new_node.start_angle_ = N_[parent_index].C_[child_index].start_angle_;
      new_node.end_angle_ = N_[parent_index].C_[child_index].end_angle_;

      std::string child_error;
      const auto child_status = new_node.generateChild(theMap.get(), stop_token, &child_error);
      if (child_status == OperationStatus::stopped || stop_token.poll())
        return stop_for_request();
      if (child_status == OperationStatus::failure)
        return fail_planning("Scoped child generation failed: " + child_error);
      if (N_.size() >= limits.max_nodes)
        return stop_for_limit(PlanningLimitReached::max_nodes);
      pending_node.emplace(std::move(new_node));
    }

    if (!pending_node)
      return fail_planning("Planner failed to construct the selected node");

    for (auto& child : pending_node->C_) {
      if (stop_token.poll())
        return stop_for_request();
      child.c_hcost_ = std::hypot(child.c_endpoint_.position.first - goal_gx,
                                  child.c_endpoint_.position.second - goal_gy);
    }
    if (stop_token.poll())
      return stop_for_request();

    N_.emplace_back(std::move(*pending_node));
    for (const auto& child : N_.back().C_) {
      if (stop_token.poll())
        return stop_for_request();
      queue.emplace_back(Candidate(child.Nindex_, child.Cindex_, child.c_gcost_ + child.c_hcost_));
      std::push_heap(queue.begin(), queue.end(), comp);
    }

    if (stop_token.poll())
      return stop_for_request();
    const auto point_location_mode = visibilityBoundaryModeForNode(*theMap, N_.back());
    const exact_geometry::Point exact_node_source =
      N_.back().is_continuous_root_
        ? exact_geometry::Point(start_gx, start_gy)
        : exact_geometry::Point(N_.back().seed_.first, N_.back().seed_.second);
    std::pair<bool, bool> b = {false, false};
    const auto point_location_status =
      detail::classifyPointInVisibilityBoundary(N_.back().visibility_region_,
                                                exact_geometry::Point(goal_gx, goal_gy),
                                                exact_node_source,
                                                point_location_mode,
                                                b,
                                                stop_token);
    if (point_location_status == OperationStatus::stopped || stop_token.poll())
      return stop_for_request();
    if (point_location_status == OperationStatus::failure)
      return fail_planning("Goal point-location validation failed");
    if (b.first || b.second) {
      bool goal_segment_crosses = false;
      if (!allow_self_crossing) {
        const auto crossing_status = hybridPathSelfCrosses({start_gx, start_gy},
                                                           N_.back().local_shortest_path_,
                                                           {goal_gx, goal_gy},
                                                           goal_segment_crosses,
                                                           stop_token);
        if (crossing_status == OperationStatus::stopped || stop_token.poll())
          return stop_for_request();
        if (crossing_status == OperationStatus::failure)
          return fail_planning("Goal path self-crossing validation failed");
      }

      if (!goal_segment_crosses) {
        ConservativeBinary64PathLength path_certificate(maximum_length_refinement_precision_);
        if (!buildPolylineLengthCertificate({start_gx, start_gy},
                                            N_.back().local_shortest_path_,
                                            {goal_gx, goal_gy},
                                            path_certificate,
                                            &stop_token)) {
          if (stop_token.poll())
            return stop_for_request();
          return fail_planning("Could not construct a path length certificate");
        }
        double certified_path_cost = std::numeric_limits<double>::quiet_NaN();
        const auto certificate_stop = [&]() { return stop_token.poll(); };
        if (!path_certificate.upperBound(certified_path_cost, certificate_stop)) {
          if (stop_token.poll())
            return stop_for_request();
          return fail_planning("Could not compute the path length upper certificate");
        }

        const size_t turning_point_count = N_.back().local_shortest_path_.size();
        // PathSolution always contains the two continuous endpoints in
        // addition to its integer turning points.  Check both the per-path
        // addition and the response-wide total before constructing/copying a
        // solution, so an adversarial search cannot grow an unbounded result.
        if (turning_point_count > limits.max_path_points - 2 ||
            accumulated_path_points > limits.max_path_points - (turning_point_count + 2)) {
          return stop_for_limit(PlanningLimitReached::max_path_points);
        }
        PathSolution solution({start_gx, start_gy},
                              N_.back().local_shortest_path_,
                              {goal_gx, goal_gy},
                              certified_path_cost,
                              N_.back().path_node_index_);
        if (stop_token.poll())
          return stop_for_request();
        result.path_solutions.emplace_back(std::move(solution));
        accumulated_path_points += turning_point_count + 2;
        if (stop_token.poll())
          return stop_for_request();
        if (objective.mode == SearchMode::top_k &&
            static_cast<int>(result.path_solutions.size()) >= objective.k) {
          result.completion = PlanningCompletion::requested_k_reached;
          break;
        }
      }
    }
  }

  if (stop_token.poll())
    return stop_for_request();

  if (queue.empty() && result.completion == PlanningCompletion::none)
    result.completion = PlanningCompletion::frontier_exhausted;

  const auto planner_end_time = PlanningClock::now();
  result.plan_time_ms =
    std::chrono::duration_cast<std::chrono::microseconds>(planner_end_time - planner_start_time)
      .count() /
    1000.0;

  result.success = !result.path_solutions.empty();
  result.expanded_nodes = N_.size();
  if (!result.success) {
    result.outcome = PlanningOutcome::no_path;
    result.message = "Planning completed but no path found";
  } else {
    result.outcome = PlanningOutcome::complete;
  }
  result.polymap = theMap;
  return result;
} catch (const std::exception& exception) {
  resetSearchState();
  PlanResult result;
  result.message = boundedPlanningExceptionMessage(exception.what());
  return result;
} catch (...) {
  resetSearchState();
  PlanResult result;
  result.message = "Planning failed with an unknown exception";
  return result;
}


}  // namespace raystar
