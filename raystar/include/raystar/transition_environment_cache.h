#pragma once

#include <raystar/polymap.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <optional>
#include <tuple>
#include <utility>
#include <vector>

namespace raystar {

// Polymap construction is not map-only: contour extraction selects the free
// component containing the root, and obstacle simplification protects the exact
// continuous root/goal positions.  Keep those inputs in the cache identity so
// a completed triangulation is reused only for an equivalent construction.
struct TransitionEnvironmentEndpoint {
  int cell_x = 0;
  int cell_y = 0;
  Point2d position = {0.0, 0.0};

  TransitionEnvironmentEndpoint() = default;

  explicit TransitionEnvironmentEndpoint(const PolymapEndpoint& endpoint)
    : cell_x(endpoint.cell_x), cell_y(endpoint.cell_y), position(endpoint.position) {}

  TransitionEnvironmentEndpoint(int x, int y, Point2d continuous_position)
    : cell_x(x), cell_y(y), position(std::move(continuous_position)) {}

  [[nodiscard]] bool operator==(const TransitionEnvironmentEndpoint& other) const noexcept {
    return cell_x == other.cell_x && cell_y == other.cell_y && position == other.position;
  }

  [[nodiscard]] bool operator<(const TransitionEnvironmentEndpoint& other) const noexcept {
    return std::tie(cell_x, cell_y, position.first, position.second) <
           std::tie(other.cell_x, other.cell_y, other.position.first, other.position.second);
  }
};

// These values either change the binary occupancy presented to Polymap or
// govern admission of its map working set.  Admission is still repeated on a
// cache hit; retaining the limits here also makes runtime parameter changes an
// explicit miss instead of silently reusing an object admitted under an older
// configuration snapshot.
struct TransitionEnvironmentPolicy {
  bool allow_unknown = false;
  int occupied_threshold = 0;
  size_t max_map_cells = 0;
  size_t max_map_bytes = 0;

  [[nodiscard]] bool operator==(const TransitionEnvironmentPolicy& other) const noexcept {
    return allow_unknown == other.allow_unknown && occupied_threshold == other.occupied_threshold &&
           max_map_cells == other.max_map_cells && max_map_bytes == other.max_map_bytes;
  }
};

class TransitionEnvironmentKey {
public:
  TransitionEnvironmentKey(std::uint64_t map_generation,
                           TransitionEnvironmentPolicy policy,
                           TransitionEnvironmentEndpoint root,
                           std::vector<TransitionEnvironmentEndpoint> goals)
    : map_generation_(map_generation)
    , policy_(policy)
    , root_(std::move(root))
    , goals_(std::move(goals)) {
    // Goal order and duplicate goal entries do not affect either reachability
    // or the protected-point membership checks used by simplification.
    std::sort(goals_.begin(), goals_.end());
    goals_.erase(std::unique(goals_.begin(), goals_.end()), goals_.end());
  }

  [[nodiscard]] bool operator==(const TransitionEnvironmentKey& other) const noexcept {
    return map_generation_ == other.map_generation_ && policy_ == other.policy_ &&
           root_ == other.root_ && goals_ == other.goals_;
  }

private:
  std::uint64_t map_generation_ = 0;
  TransitionEnvironmentPolicy policy_;
  TransitionEnvironmentEndpoint root_;
  std::vector<TransitionEnvironmentEndpoint> goals_;
};

// A capacity-one strong cache is sufficient for the TMV GCP -> UPS pipeline
// and gives a hard memory bound.  Shared ownership keeps an environment alive
// for an in-flight UPS request even when a map callback invalidates the entry.
// Only completed, const Polymaps are accepted; no Core search-tree state or
// request StopToken is retained.
class CompletedTransitionEnvironmentCache {
public:
  [[nodiscard]] std::shared_ptr<const Polymap> find(const TransitionEnvironmentKey& key) const {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!entry_ || !(entry_->key == key))
      return {};
    return entry_->environment;
  }

  void store(TransitionEnvironmentKey key, std::shared_ptr<const Polymap> environment) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (!environment) {
      entry_.reset();
      return;
    }
    entry_.emplace(Entry{std::move(key), std::move(environment)});
  }

  void clear() {
    std::lock_guard<std::mutex> lock(mutex_);
    entry_.reset();
  }

  [[nodiscard]] bool empty() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return !entry_.has_value();
  }

private:
  struct Entry {
    TransitionEnvironmentKey key;
    std::shared_ptr<const Polymap> environment;
  };

  mutable std::mutex mutex_;
  std::optional<Entry> entry_;
};

}  // namespace raystar
