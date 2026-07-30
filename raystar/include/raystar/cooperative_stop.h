#pragma once

#include <functional>
#include <utility>

namespace raystar {

using StopPredicate = std::function<bool()>;

enum class OperationStatus { success, failure, stopped };

// A lightweight cooperative-stop token.  poll() evaluates the predicate until
// it first requests a stop, then latches that state so all later polls remain
// stopped without invoking the predicate again.
class StopToken {
public:
  StopToken() = default;

  explicit StopToken(StopPredicate predicate) : predicate_(std::move(predicate)) {}

  [[nodiscard]] bool poll() const {
    if (!requested_ && predicate_ && predicate_())
      requested_ = true;
    return requested_;
  }

  [[nodiscard]] bool requested() const noexcept {
    return requested_;
  }

private:
  StopPredicate predicate_;
  mutable bool requested_ = false;
};

}  // namespace raystar
