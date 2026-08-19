#pragma once

#include <chrono>
#include <cstddef>
#include <limits>
#include <string>

#include <raystar/cooperative_stop.h>

namespace raystar {

// Request-wide limits shared by the Core, ROS boundary, and Polymap factory.
// Keeping this type independent of raystar_core.h prevents Polymap's
// resource-admission API from depending on the search-tree declarations.
struct PlanningLimits {
  static constexpr size_t kDefaultMaxMapCells = 8u * 1024u * 1024u;
  static constexpr size_t kDefaultMaxMapBytes = 512u * 1024u * 1024u;
  static constexpr size_t kDefaultMaxPathPoints = 100000u;
  // Debug-tree serialization is opt-in.  Keeping the default at zero avoids
  // retaining or exporting a potentially large tree for normal plans.
  static constexpr size_t kDefaultMaxDebugNodes = 0u;
  static constexpr size_t kDefaultMaxResponseBytes = 64u * 1024u * 1024u;
  static constexpr size_t kDefaultMaxCostBoundedPaths = 1000u;
  static constexpr size_t kDefaultMaxMultiGoalCount = 128u;
  static constexpr size_t kDefaultMaxTransitionConfigurations = 4096u;
  static constexpr size_t kDefaultMaxTransitionPairs = 1000u;

  int max_k = std::numeric_limits<int>::max();
  // Hard safety cap for exhaustive cost-bounded enumeration. Reaching this
  // cap is a partial-search limit unless the frontier has already crossed the
  // requested cost bound.
  size_t max_cost_bounded_paths = kDefaultMaxCostBoundedPaths;
  // Admission bound for one shared-tree multi-goal request. The per-goal
  // path cap above still applies independently, while max_path_points bounds
  // aggregate retained path geometry across the complete request.
  size_t max_multi_goal_count = kDefaultMaxMultiGoalCount;
  // Independent UPS batch admission limits. Configuration count is the
  // flattened union of all layers; pair count is the explicit directed edge
  // array. They are intentionally separate from GCP goal/path limits.
  size_t max_transition_configurations = kDefaultMaxTransitionConfigurations;
  size_t max_transition_pairs = kDefaultMaxTransitionPairs;
  // Counts fully constructed Nodes, including the root Node.
  size_t max_nodes = static_cast<size_t>(std::numeric_limits<int>::max());
  // Cooperative deadline. Zero requests an immediate timeout; max() disables
  // the deadline for direct Core callers.
  std::chrono::milliseconds planning_timeout = std::chrono::milliseconds::max();
  // Optional cooperative-cancellation hook. Callers that signal it from
  // another thread must provide a thread-safe predicate.
  StopPredicate cancel_requested;
  // Fixed map admission is checked before either Core or Polymap copies
  // GridMap. Reference-shortening Polymaps additionally charge unsimplified
  // contour topology against the unused byte budget before constructing CDT.
  size_t max_map_cells = kDefaultMaxMapCells;
  size_t max_map_bytes = kDefaultMaxMapBytes;
  // Aggregate structural/output limits used by Core and the ROS adapter.
  size_t max_path_points = kDefaultMaxPathPoints;
  size_t max_debug_nodes = kDefaultMaxDebugNodes;
  size_t max_response_bytes = kDefaultMaxResponseBytes;
};

struct MapResourceEstimate {
  size_t cell_count = 0;
  size_t estimated_bytes = 0;
};

// Fixed whole-map storage that can coexist during construction includes the
// Core work map, Polymap copy, flood-fill mask, and old/new vertex registries.
inline constexpr size_t kEstimatedPlannerMapBytesPerCell = 32u;

/// Validate map shape and its resource admission budget without allocating.
[[nodiscard]] bool validateMapResourceBudget(size_t width,
                                             size_t height,
                                             size_t data_size,
                                             const PlanningLimits& limits,
                                             MapResourceEstimate& estimate,
                                             std::string& error);

}  // namespace raystar
