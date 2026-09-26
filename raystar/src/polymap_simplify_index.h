#pragma once

// Candidate index for the polygon-simplification safety checks
// (SIMPLIFY_STACK_SCAN_DESIGN.md).  A uniform bucket grid over the obstacle
// bounding box: vertices are registered in the bucket containing them and
// every edge is registered in every bucket its axis-aligned bounding box
// overlaps.  Queries return every registration inside the buckets overlapped
// by the query box.
//
// Completeness lemma (design doc, Lemma 1): any vertex/edge that is
// non-disjoint from a query object shares a point p with it; p lies in the
// edge's bbox and in the query bbox, so the bucket containing p is both
// registered and queried.  Pure integer interval arithmetic -- no geometric
// predicates are involved, so the candidate set is a conservative superset
// by construction and the exact EPECK predicates run unchanged on it.
//
// The index identifies vertices and edges by their obstacle index and
// ORIGINAL vertex id (stable for the lifetime of one simplification call;
// the caller keeps ring connectivity in original-id space and compacts the
// ring vector only after the call).  Edge id e = the segment from original
// vertex e to the caller's current successor of e.
//
// Build refuses inputs that would change the legacy semantics or the
// bookkeeping assumptions -- a duplicate coordinate anywhere or a
// zero-length edge (the legacy chord check carries a global sentinel for
// those, polymap.cpp).  The caller must then run the full-scan candidate
// provider instead.  Production inputs never trip this: every simplify call
// is preceded by validateObstacleTopology, which rejects both; the gate
// exists for the incremental Stage 7 unfreeze inputs, which may legally
// violate the preconditions.

#include <cstdint>
#include <utility>
#include <vector>

namespace raystar {

class Obs;

namespace polymap_impl {

class SimplifyCandidateIndex {
public:
  struct VertexRecord {
    int obstacle;
    int vertex;  // original vertex id
  };
  struct EdgeRecord {
    int obstacle;
    int edge;  // original id of the edge's source vertex
  };

  // bucket_size is in grid units per bucket side.
  explicit SimplifyCandidateIndex(int bucket_size = 8);

  // Registers every vertex and every edge of every ring (empty rings --
  // tombstones -- contribute nothing).  Returns false and leaves the index
  // unusable when the input is not indexable (duplicate coordinate or
  // zero-length edge anywhere).
  [[nodiscard]] bool build(const std::vector<Obs>& obstacles);

  // Append every registration whose bucket range overlaps the closed box
  // [min_x, max_x] x [min_y, max_y].  Each record is reported once per call
  // (epoch-based dedup).  The box may extend beyond the indexed area.
  void verticesInBox(
    int min_x, int min_y, int max_x, int max_y, std::vector<VertexRecord>& out) const;
  void edgesInBox(int min_x, int min_y, int max_x, int max_y, std::vector<EdgeRecord>& out) const;

  // Maintenance for one vertex removal (obstacle o, original ids
  // previous/current/next with their CURRENT geometry): drops the removed
  // vertex registration, drops the two replaced edge registrations
  // (previous->current, current->next) and registers the chord
  // (previous->next) under edge id `previous`.
  void applyRemoval(int obstacle,
                    int previous,
                    const std::pair<int, int>& previous_xy,
                    int current,
                    const std::pair<int, int>& current_xy,
                    int next,
                    const std::pair<int, int>& next_xy);

private:
  struct Bucket {
    std::vector<VertexRecord> vertices;
    std::vector<EdgeRecord> edges;
  };

  [[nodiscard]] int bucketOf(int x, int y) const;
  void bucketRange(int min_x,
                   int min_y,
                   int max_x,
                   int max_y,
                   int& bucket_min_x,
                   int& bucket_min_y,
                   int& bucket_max_x,
                   int& bucket_max_y) const;
  void registerEdge(const EdgeRecord& record,
                    const std::pair<int, int>& from,
                    const std::pair<int, int>& to);
  void unregisterEdge(const EdgeRecord& record,
                      const std::pair<int, int>& from,
                      const std::pair<int, int>& to);

  int bucket_size_;
  int origin_x_ = 0;
  int origin_y_ = 0;
  int buckets_x_ = 0;
  int buckets_y_ = 0;
  std::vector<Bucket> buckets_;
  // Per-record dedup epochs, addressed by (obstacle, original id).
  mutable std::vector<std::vector<std::uint32_t>> vertex_epoch_;
  mutable std::vector<std::vector<std::uint32_t>> edge_epoch_;
  mutable std::uint32_t epoch_ = 0;
};

}  // namespace polymap_impl
}  // namespace raystar
