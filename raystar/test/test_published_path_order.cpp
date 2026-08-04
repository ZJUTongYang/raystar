#include "published_path_order.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace raystar {
namespace {

TEST(PublishedPathOrder, SortsByFinalPublicCostAndPreservesStableSourceAssociation) {
  std::vector<raystar_interfaces::msg::PathResult> paths(3);
  // Model a legal Core order whose first two records reverse after independent
  // serialized-world length ceilings are finalized.
  paths[0].cost = std::nextafter(1.0, std::numeric_limits<double>::infinity());
  paths[0].topology_path.header.frame_id = "higher-first";
  paths[1].cost = 1.0;
  paths[1].topology_path.header.frame_id = "lower-second";
  paths[2].cost = 1.0;
  paths[2].topology_path.header.frame_id = "equal-third";
  std::vector<std::size_t> source_indices{10u, 20u, 30u};

  ASSERT_TRUE(stableSortPublishedPathsWithSources(paths, source_indices));

  ASSERT_EQ(paths.size(), 3u);
  EXPECT_DOUBLE_EQ(paths[0].cost, 1.0);
  EXPECT_DOUBLE_EQ(paths[1].cost, 1.0);
  EXPECT_GT(paths[2].cost, paths[1].cost);
  EXPECT_EQ(paths[0].topology_path.header.frame_id, "lower-second");
  EXPECT_EQ(paths[1].topology_path.header.frame_id, "equal-third");
  EXPECT_EQ(paths[2].topology_path.header.frame_id, "higher-first");
  EXPECT_EQ(source_indices, (std::vector<std::size_t>{20u, 30u, 10u}));
}

TEST(PublishedPathOrder, RejectsMismatchedParallelSourceCardinalityWithoutMutation) {
  std::vector<raystar_interfaces::msg::PathResult> paths(1);
  paths.front().cost = 1.0;
  std::vector<std::size_t> source_indices;

  EXPECT_FALSE(stableSortPublishedPathsWithSources(paths, source_indices));
  ASSERT_EQ(paths.size(), 1u);
  EXPECT_DOUBLE_EQ(paths.front().cost, 1.0);
  EXPECT_TRUE(source_indices.empty());
}

}  // namespace
}  // namespace raystar
