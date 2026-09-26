#include "polymap_simplify_index.h"

#include <raystar/polymap.h>

#include <algorithm>
#include <unordered_set>

namespace raystar {
namespace polymap_impl {

namespace {

// Packs an integer coordinate pair for the duplicate-coordinate gate.
inline std::uint64_t packCoordinate(int x, int y) {
  return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(x)) << 32) |
         static_cast<std::uint64_t>(static_cast<std::uint32_t>(y));
}

}  // namespace

SimplifyCandidateIndex::SimplifyCandidateIndex(int bucket_size)
  : bucket_size_(bucket_size < 1 ? 1 : bucket_size) {}

int SimplifyCandidateIndex::bucketOf(int x, int y) const {
  int bucket_x = (x - origin_x_) / bucket_size_;
  int bucket_y = (y - origin_y_) / bucket_size_;
  bucket_x = std::clamp(bucket_x, 0, buckets_x_ - 1);
  bucket_y = std::clamp(bucket_y, 0, buckets_y_ - 1);
  return bucket_y * buckets_x_ + bucket_x;
}

void SimplifyCandidateIndex::bucketRange(int min_x,
                                         int min_y,
                                         int max_x,
                                         int max_y,
                                         int& bucket_min_x,
                                         int& bucket_min_y,
                                         int& bucket_max_x,
                                         int& bucket_max_y) const {
  bucket_min_x = std::clamp((min_x - origin_x_) / bucket_size_, 0, buckets_x_ - 1);
  bucket_min_y = std::clamp((min_y - origin_y_) / bucket_size_, 0, buckets_y_ - 1);
  bucket_max_x = std::clamp((max_x - origin_x_) / bucket_size_, 0, buckets_x_ - 1);
  bucket_max_y = std::clamp((max_y - origin_y_) / bucket_size_, 0, buckets_y_ - 1);
}

bool SimplifyCandidateIndex::build(const std::vector<Obs>& obstacles) {
  buckets_.clear();
  vertex_epoch_.clear();
  edge_epoch_.clear();
  epoch_ = 0;

  // Pass 1: bounding box, the indexability gate, and the epoch tables.
  bool any_vertex = false;
  int min_x = 0;
  int min_y = 0;
  int max_x = 0;
  int max_y = 0;
  std::unordered_set<std::uint64_t> seen_coordinates;
  vertex_epoch_.resize(obstacles.size());
  edge_epoch_.resize(obstacles.size());
  for (size_t obstacle = 0; obstacle < obstacles.size(); ++obstacle) {
    const auto& ring = obstacles[obstacle].ordered_vertices_;
    vertex_epoch_[obstacle].assign(ring.size(), 0);
    edge_epoch_[obstacle].assign(ring.size(), 0);
    for (size_t vertex = 0; vertex < ring.size(); ++vertex) {
      const auto& coordinate = ring[vertex];
      const auto& next = ring[(vertex + 1) % ring.size()];
      if (coordinate == next)
        return false;  // zero-length edge: legacy sentinel semantics apply
      if (!seen_coordinates.insert(packCoordinate(coordinate.first, coordinate.second)).second)
        return false;  // duplicate coordinate anywhere (Stage 7 unfreeze inputs)
      if (!any_vertex) {
        min_x = max_x = coordinate.first;
        min_y = max_y = coordinate.second;
        any_vertex = true;
      } else {
        min_x = std::min(min_x, coordinate.first);
        max_x = std::max(max_x, coordinate.first);
        min_y = std::min(min_y, coordinate.second);
        max_y = std::max(max_y, coordinate.second);
      }
    }
  }
  if (!any_vertex)
    return false;  // nothing to index; callers with no rings never query

  origin_x_ = min_x;
  origin_y_ = min_y;
  buckets_x_ = (max_x - min_x) / bucket_size_ + 1;
  buckets_y_ = (max_y - min_y) / bucket_size_ + 1;
  buckets_.assign(static_cast<size_t>(buckets_x_) * static_cast<size_t>(buckets_y_), Bucket{});

  // Pass 2: registrations.
  for (size_t obstacle = 0; obstacle < obstacles.size(); ++obstacle) {
    const auto& ring = obstacles[obstacle].ordered_vertices_;
    for (size_t vertex = 0; vertex < ring.size(); ++vertex) {
      const auto& coordinate = ring[vertex];
      buckets_[static_cast<size_t>(bucketOf(coordinate.first, coordinate.second))]
        .vertices.push_back(VertexRecord{static_cast<int>(obstacle), static_cast<int>(vertex)});
      registerEdge(EdgeRecord{static_cast<int>(obstacle), static_cast<int>(vertex)},
                   coordinate,
                   ring[(vertex + 1) % ring.size()]);
    }
  }
  return true;
}

void SimplifyCandidateIndex::registerEdge(const EdgeRecord& record,
                                          const std::pair<int, int>& from,
                                          const std::pair<int, int>& to) {
  int bucket_min_x, bucket_min_y, bucket_max_x, bucket_max_y;
  bucketRange(std::min(from.first, to.first),
              std::min(from.second, to.second),
              std::max(from.first, to.first),
              std::max(from.second, to.second),
              bucket_min_x,
              bucket_min_y,
              bucket_max_x,
              bucket_max_y);
  for (int bucket_y = bucket_min_y; bucket_y <= bucket_max_y; ++bucket_y)
    for (int bucket_x = bucket_min_x; bucket_x <= bucket_max_x; ++bucket_x)
      buckets_[static_cast<size_t>(bucket_y) * static_cast<size_t>(buckets_x_) +
               static_cast<size_t>(bucket_x)]
        .edges.push_back(record);
}

void SimplifyCandidateIndex::unregisterEdge(const EdgeRecord& record,
                                            const std::pair<int, int>& from,
                                            const std::pair<int, int>& to) {
  int bucket_min_x, bucket_min_y, bucket_max_x, bucket_max_y;
  bucketRange(std::min(from.first, to.first),
              std::min(from.second, to.second),
              std::max(from.first, to.first),
              std::max(from.second, to.second),
              bucket_min_x,
              bucket_min_y,
              bucket_max_x,
              bucket_max_y);
  for (int bucket_y = bucket_min_y; bucket_y <= bucket_max_y; ++bucket_y)
    for (int bucket_x = bucket_min_x; bucket_x <= bucket_max_x; ++bucket_x) {
      auto& edges = buckets_[static_cast<size_t>(bucket_y) * static_cast<size_t>(buckets_x_) +
                             static_cast<size_t>(bucket_x)]
                      .edges;
      for (size_t index = 0; index < edges.size(); ++index) {
        if (edges[index].obstacle == record.obstacle && edges[index].edge == record.edge) {
          edges[index] = edges.back();
          edges.pop_back();
          break;
        }
      }
    }
}

void SimplifyCandidateIndex::advanceEpoch() const {
  ++epoch_;
  if (epoch_ == 0) {
    // uint32 wrap after ~4e9 queries: a record stamped before the wrap could
    // alias the fresh counter and be silently deduped away -- a missed
    // candidate.  Re-zero the tables and restart at 1.
    for (auto& table : vertex_epoch_) std::fill(table.begin(), table.end(), 0);
    for (auto& table : edge_epoch_) std::fill(table.begin(), table.end(), 0);
    epoch_ = 1;
  }
}

void SimplifyCandidateIndex::verticesInBox(
  int min_x, int min_y, int max_x, int max_y, std::vector<VertexRecord>& out) const {
  advanceEpoch();
  int bucket_min_x, bucket_min_y, bucket_max_x, bucket_max_y;
  bucketRange(min_x, min_y, max_x, max_y, bucket_min_x, bucket_min_y, bucket_max_x, bucket_max_y);
  for (int bucket_y = bucket_min_y; bucket_y <= bucket_max_y; ++bucket_y)
    for (int bucket_x = bucket_min_x; bucket_x <= bucket_max_x; ++bucket_x) {
      const auto& bucket =
        buckets_[static_cast<size_t>(bucket_y) * static_cast<size_t>(buckets_x_) +
                 static_cast<size_t>(bucket_x)];
      for (const auto& record : bucket.vertices) {
        auto& seen =
          vertex_epoch_[static_cast<size_t>(record.obstacle)][static_cast<size_t>(record.vertex)];
        if (seen == epoch_)
          continue;
        seen = epoch_;
        out.push_back(record);
      }
    }
}

void SimplifyCandidateIndex::edgesInBox(
  int min_x, int min_y, int max_x, int max_y, std::vector<EdgeRecord>& out) const {
  advanceEpoch();
  int bucket_min_x, bucket_min_y, bucket_max_x, bucket_max_y;
  bucketRange(min_x, min_y, max_x, max_y, bucket_min_x, bucket_min_y, bucket_max_x, bucket_max_y);
  for (int bucket_y = bucket_min_y; bucket_y <= bucket_max_y; ++bucket_y)
    for (int bucket_x = bucket_min_x; bucket_x <= bucket_max_x; ++bucket_x) {
      const auto& bucket =
        buckets_[static_cast<size_t>(bucket_y) * static_cast<size_t>(buckets_x_) +
                 static_cast<size_t>(bucket_x)];
      for (const auto& record : bucket.edges) {
        auto& seen =
          edge_epoch_[static_cast<size_t>(record.obstacle)][static_cast<size_t>(record.edge)];
        if (seen == epoch_)
          continue;
        seen = epoch_;
        out.push_back(record);
      }
    }
}

void SimplifyCandidateIndex::applyRemoval(int obstacle,
                                          int previous,
                                          const std::pair<int, int>& previous_xy,
                                          int current,
                                          const std::pair<int, int>& current_xy,
                                          int next,
                                          const std::pair<int, int>& next_xy) {
  (void)next;  // identifies the successor for the caller; the chord keeps `previous` as id
  // Drop the removed vertex registration.
  auto& vertices =
    buckets_[static_cast<size_t>(bucketOf(current_xy.first, current_xy.second))].vertices;
  for (size_t index = 0; index < vertices.size(); ++index) {
    if (vertices[index].obstacle == obstacle && vertices[index].vertex == current) {
      vertices[index] = vertices.back();
      vertices.pop_back();
      break;
    }
  }
  // Replace the two edges with the chord, registered under the surviving
  // source vertex id.
  unregisterEdge(EdgeRecord{obstacle, previous}, previous_xy, current_xy);
  unregisterEdge(EdgeRecord{obstacle, current}, current_xy, next_xy);
  registerEdge(EdgeRecord{obstacle, previous}, previous_xy, next_xy);
}

}  // namespace polymap_impl
}  // namespace raystar
