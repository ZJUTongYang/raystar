#include <gtest/gtest.h>

#include <raystar/transition_environment_cache.h>

#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace raystar {
namespace {

std::shared_ptr<const Polymap> makeCompletedEnvironment() {
  GridMap map;
  map.width = 8;
  map.height = 8;
  map.data.assign(static_cast<size_t>(map.width) * map.height, 0);
  for (unsigned int x = 0; x < map.width; ++x) {
    map(x, 0) = 1;
    map(x, map.height - 1) = 1;
  }
  for (unsigned int y = 0; y < map.height; ++y) {
    map(0, y) = 1;
    map(map.width - 1, y) = 1;
  }

  auto created = Polymap::create(map, 2, 2, 5, 5, Point2d{2.5, 2.5}, Point2d{5.5, 5.5});
  if (!created) {
    throw std::runtime_error("Test Polymap construction failed: " + created.error);
  }
  return std::make_shared<Polymap>(std::move(*created.value));
}

TransitionEnvironmentPolicy makePolicy() {
  return TransitionEnvironmentPolicy{false, 99, 1'000'000, 64'000'000};
}

TransitionEnvironmentEndpoint makeEndpoint(int cell_x, int cell_y, double x, double y) {
  return TransitionEnvironmentEndpoint(cell_x, cell_y, Point2d{x, y});
}

TEST(CompletedTransitionEnvironmentCacheTest, CanonicalizesGoalSetButMatchesEveryOtherInput) {
  const auto environment = makeCompletedEnvironment();
  CompletedTransitionEnvironmentCache cache;
  const auto base = makeEndpoint(2, 2, 2.5, 2.5);
  const auto first_goal = makeEndpoint(5, 5, 5.5, 5.5);
  const auto second_goal = makeEndpoint(4, 5, 4.5, 5.5);
  const TransitionEnvironmentKey stored_key(
    7, makePolicy(), base, {second_goal, first_goal, first_goal});
  cache.store(stored_key, environment);

  const TransitionEnvironmentKey reordered_key(7, makePolicy(), base, {first_goal, second_goal});
  EXPECT_EQ(cache.find(reordered_key), environment);

  auto different_policy = makePolicy();
  different_policy.allow_unknown = true;
  EXPECT_FALSE(
    cache.find(TransitionEnvironmentKey(7, different_policy, base, {first_goal, second_goal})));
  different_policy = makePolicy();
  different_policy.occupied_threshold = 50;
  EXPECT_FALSE(
    cache.find(TransitionEnvironmentKey(7, different_policy, base, {first_goal, second_goal})));
  different_policy = makePolicy();
  --different_policy.max_map_cells;
  EXPECT_FALSE(
    cache.find(TransitionEnvironmentKey(7, different_policy, base, {first_goal, second_goal})));
  different_policy = makePolicy();
  --different_policy.max_map_bytes;
  EXPECT_FALSE(
    cache.find(TransitionEnvironmentKey(7, different_policy, base, {first_goal, second_goal})));

  EXPECT_FALSE(
    cache.find(TransitionEnvironmentKey(8, makePolicy(), base, {first_goal, second_goal})));
  EXPECT_FALSE(cache.find(TransitionEnvironmentKey(
    7, makePolicy(), makeEndpoint(2, 2, 2.25, 2.5), {first_goal, second_goal})));
  EXPECT_FALSE(cache.find(
    TransitionEnvironmentKey(7, makePolicy(), base, {makeEndpoint(5, 5, 5.25, 5.5), second_goal})));
  EXPECT_FALSE(cache.find(TransitionEnvironmentKey(7, makePolicy(), base, {first_goal})));
}

TEST(CompletedTransitionEnvironmentCacheTest, ReplacementIsCapacityOneAndOwnershipIsRaii) {
  CompletedTransitionEnvironmentCache cache;
  const auto base = makeEndpoint(2, 2, 2.5, 2.5);
  const auto goal = makeEndpoint(5, 5, 5.5, 5.5);
  const TransitionEnvironmentKey first_key(1, makePolicy(), base, {goal});
  const TransitionEnvironmentKey second_key(2, makePolicy(), base, {goal});

  auto first = makeCompletedEnvironment();
  std::weak_ptr<const Polymap> first_weak = first;
  cache.store(first_key, first);
  first.reset();
  EXPECT_FALSE(first_weak.expired());

  auto second = makeCompletedEnvironment();
  std::weak_ptr<const Polymap> second_weak = second;
  cache.store(second_key, second);
  EXPECT_TRUE(first_weak.expired());
  EXPECT_FALSE(cache.find(first_key));

  auto in_flight = cache.find(second_key);
  ASSERT_EQ(in_flight, second);
  second.reset();
  cache.clear();
  EXPECT_TRUE(cache.empty());
  EXPECT_FALSE(second_weak.expired());
  in_flight.reset();
  EXPECT_TRUE(second_weak.expired());
}

}  // namespace
}  // namespace raystar
