#pragma once

#include <raystar_interfaces/msg/path_result.hpp>

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace raystar {

// Core ranks candidates by grid-space cost, while the public cost is the
// conservative ceiling of the two serialized world-coordinate polylines.
// Affine-transform and interpolation rounding can change that last ULP, so
// restore the public nondecreasing-cost contract after finalization.  Move the
// parallel source identifiers through the same stable permutation to preserve
// the response/visualization association and Core tie order.
template<typename SourceIdentifier>
bool stableSortPublishedPathsWithSources(
  std::vector<raystar_interfaces::msg::PathResult>& paths,
  std::vector<SourceIdentifier>& sources) {
  if (paths.size() != sources.size())
    return false;

  std::vector<std::size_t> order;
  order.reserve(paths.size());
  for (std::size_t index = 0; index < paths.size(); ++index)
    order.emplace_back(index);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
    return paths[lhs].cost < paths[rhs].cost;
  });

  std::vector<raystar_interfaces::msg::PathResult> sorted_paths;
  std::vector<SourceIdentifier> sorted_sources;
  sorted_paths.reserve(paths.size());
  sorted_sources.reserve(sources.size());
  for (const std::size_t index : order) {
    sorted_paths.emplace_back(std::move(paths[index]));
    sorted_sources.emplace_back(std::move(sources[index]));
  }
  paths = std::move(sorted_paths);
  sources = std::move(sorted_sources);
  return true;
}

}  // namespace raystar
