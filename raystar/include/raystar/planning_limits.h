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

  int max_k = std::numeric_limits<int>::max();
  // Counts fully constructed Nodes, including the root Node.
  size_t max_nodes = static_cast<size_t>(std::numeric_limits<int>::max());
  // Cooperative deadline. Zero requests an immediate timeout; max() disables
  // the deadline for direct Core callers.
  std::chrono::milliseconds planning_timeout = std::chrono::milliseconds::max();
  // Optional cooperative-cancellation hook. Callers that signal it from
  // another thread must provide a thread-safe predicate.
  StopPredicate cancel_requested;
  // Map admission is checked before either Core or Polymap copies GridMap.
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
