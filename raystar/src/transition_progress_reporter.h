#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>

namespace raystar {

// The ROS adapter owns the callback and keeps it alive for the complete
// synchronous executeTransitionPlanning() call. The callback itself may
// allocate or throw; TransitionProgressReporter contains that failure so
// feedback can never change the planning result or unwind the Action worker.
using TransitionProgressCallback =
  std::function<void(std::uint32_t requested, std::uint32_t completed, const char* stage)>;

// Bound progress traffic for large TMV transition graphs while preserving
// useful progress across the whole batch. A normal completed request emits
// four stage records, at most 95 evenly distributed pair records (including
// pair 1 and pair N), and one terminal record: never more than 100 messages.
class TransitionProgressReporter {
public:
  static constexpr std::size_t kMaxFeedbackMessages = 100;
  static constexpr std::size_t kMaxPairFeedbackMessages = 95;

  TransitionProgressReporter(std::size_t requested,
                             const TransitionProgressCallback& callback) noexcept
    : requested_(boundedUint32(requested)),
      pair_feedback_count_(
        std::min<std::size_t>(requested_, kMaxPairFeedbackMessages)),
      callback_(&callback) {}

  void publishStage(const char* stage) noexcept {
    emit(last_observed_completed_, stage, false);
  }

  void publishPairCompletion(std::size_t completed, const char* stage) noexcept {
    if (requested_ == 0 || next_pair_feedback_ >= pair_feedback_count_)
      return;

    const std::uint32_t bounded_completed =
      std::min(requested_, boundedUint32(completed));
    if (bounded_completed < last_observed_completed_)
      return;
    last_observed_completed_ = bounded_completed;

    if (bounded_completed < pairFeedbackTarget(next_pair_feedback_))
      return;

    emit(bounded_completed, stage, false);
    do {
      ++next_pair_feedback_;
    } while (next_pair_feedback_ < pair_feedback_count_ &&
             pairFeedbackTarget(next_pair_feedback_) <= bounded_completed);
  }

  void publishFinal(const char* stage) noexcept {
    if (final_published_)
      return;
    final_published_ = true;
    last_observed_completed_ = requested_;
    emit(requested_, stage, true);
  }

private:
  static std::uint32_t boundedUint32(std::size_t value) noexcept {
    return value > std::numeric_limits<std::uint32_t>::max()
             ? std::numeric_limits<std::uint32_t>::max()
             : static_cast<std::uint32_t>(value);
  }

  std::uint32_t pairFeedbackTarget(std::size_t index) const noexcept {
    if (pair_feedback_count_ <= 1)
      return requested_;
    const std::uint64_t numerator =
      static_cast<std::uint64_t>(index) * static_cast<std::uint64_t>(requested_ - 1);
    return 1u + static_cast<std::uint32_t>(
                  numerator / static_cast<std::uint64_t>(pair_feedback_count_ - 1));
  }

  void emit(std::uint32_t completed, const char* stage, bool final) noexcept {
    if (callback_ == nullptr || !*callback_)
      return;
    // Reserve the hundredth message for publishFinal(), even if a caller adds
    // an extra diagnostic stage in the future.
    const std::size_t limit = final ? kMaxFeedbackMessages : kMaxFeedbackMessages - 1;
    if (emitted_messages_ >= limit)
      return;
    completed = std::max(completed, last_published_completed_);
    last_published_completed_ = completed;
    ++emitted_messages_;
    try {
      (*callback_)(requested_, completed, stage == nullptr ? "" : stage);
    } catch (...) {
      // Feedback is best-effort and must never affect planning/cancellation.
    }
  }

  std::uint32_t requested_ = 0;
  std::uint32_t last_observed_completed_ = 0;
  std::uint32_t last_published_completed_ = 0;
  std::size_t pair_feedback_count_ = 0;
  std::size_t next_pair_feedback_ = 0;
  std::size_t emitted_messages_ = 0;
  bool final_published_ = false;
  const TransitionProgressCallback* callback_ = nullptr;
};

}  // namespace raystar
