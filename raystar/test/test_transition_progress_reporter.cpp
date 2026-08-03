#include "transition_progress_reporter.h"

#include <gtest/gtest.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

struct ProgressRecord {
  std::uint32_t requested = 0;
  std::uint32_t completed = 0;
  std::string stage;
};

TEST(TransitionProgressReporterTest, ThrottlesLargeBatchAndPreservesEndpoints) {
  constexpr std::uint32_t requested = 649;
  std::vector<ProgressRecord> records;
  const raystar::TransitionProgressCallback callback =
    [&](std::uint32_t total, std::uint32_t completed, const char* stage) {
      records.push_back({total, completed, stage});
    };
  raystar::TransitionProgressReporter reporter(requested, callback);

  reporter.publishStage("validating transition request");
  reporter.publishStage("preparing transition environment");
  reporter.publishStage("validating tether configurations");
  reporter.publishStage("shortening transition pairs");
  for (std::uint32_t completed = 1; completed <= requested; ++completed)
    reporter.publishPairCompletion(completed, "shortening transition pairs");
  reporter.publishFinal("transition batch complete");

  ASSERT_FALSE(records.empty());
  EXPECT_EQ(records.size(), raystar::TransitionProgressReporter::kMaxFeedbackMessages);
  EXPECT_EQ(records.front().requested, requested);
  EXPECT_EQ(records.front().completed, 0u);
  EXPECT_EQ(records.front().stage, "validating transition request");
  ASSERT_GE(records.size(), 4u);
  EXPECT_EQ(records[1].completed, 0u);
  EXPECT_EQ(records[1].stage, "preparing transition environment");
  EXPECT_EQ(records[2].completed, 0u);
  EXPECT_EQ(records[2].stage, "validating tether configurations");
  EXPECT_EQ(records[3].completed, 0u);
  EXPECT_EQ(records[3].stage, "shortening transition pairs");
  EXPECT_EQ(records.back().requested, requested);
  EXPECT_EQ(records.back().completed, requested);
  EXPECT_EQ(records.back().stage, "transition batch complete");

  std::uint32_t previous_completed = 0;
  bool saw_first_pair = false;
  bool saw_last_pair = false;
  for (const auto& record : records) {
    EXPECT_EQ(record.requested, requested);
    EXPECT_GE(record.completed, previous_completed);
    EXPECT_FALSE(record.stage.empty());
    previous_completed = record.completed;
    const bool pair_record = record.stage == "shortening transition pairs";
    saw_first_pair = saw_first_pair || (pair_record && record.completed == 1u);
    saw_last_pair = saw_last_pair || (pair_record && record.completed == requested);
  }
  EXPECT_TRUE(saw_first_pair);
  EXPECT_TRUE(saw_last_pair);
}

TEST(TransitionProgressReporterTest, ContainsThrowingCallbackAndKeepsFinalProgress) {
  std::vector<ProgressRecord> records;
  std::size_t callback_count = 0;
  const raystar::TransitionProgressCallback callback =
    [&](std::uint32_t requested, std::uint32_t completed, const char* stage) {
      ++callback_count;
      if (callback_count == 2)
        throw std::runtime_error("injected feedback failure");
      records.push_back({requested, completed, stage});
    };
  raystar::TransitionProgressReporter reporter(3, callback);

  EXPECT_NO_THROW(reporter.publishStage("validating transition request"));
  EXPECT_NO_THROW(reporter.publishStage("preparing transition environment"));
  EXPECT_NO_THROW(reporter.publishStage("validating tether configurations"));
  EXPECT_NO_THROW(reporter.publishStage("shortening transition pairs"));
  for (std::uint32_t completed = 1; completed <= 3; ++completed)
    EXPECT_NO_THROW(reporter.publishPairCompletion(completed, "shortening transition pairs"));
  EXPECT_NO_THROW(reporter.publishFinal("transition batch complete"));

  ASSERT_FALSE(records.empty());
  EXPECT_EQ(records.front().completed, 0u);
  EXPECT_EQ(records.back().completed, 3u);
  EXPECT_EQ(records.back().stage, "transition batch complete");
}

}  // namespace
