#pragma once

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <vector>

namespace raystar {

// Accumulate an interval for the mathematical Euclidean length of a polyline
// whose endpoints are serialized binary64 coordinates.  Segment squares and
// every interval comparison are exact dyadic-integer operations. std::hypot
// is used only to seed a proved segment bracket; no result depends on the
// process rounding mode or on the quality of that seed.
//
// A sum of per-segment binary64 ceilings is generally not the tight ceiling
// of the mathematical polyline length. upperBound() therefore refines the
// complete sum of square roots and returns the unique smallest binary64 value
// which covers it. lowerBound() is deliberately cheaper: it sums proved
// segment floors exactly and rounds that sum downward. compareTo() adaptively
// refines the complete radical sum for strict inclusive-bound decisions.
class ConservativeBinary64PathLength {
public:
  enum class Comparison { less, equal, greater };

  explicit ConservativeBinary64PathLength(
    unsigned int maximum_refinement_precision = 16384)
    : maximum_refinement_precision_(
        std::min(maximum_refinement_precision, kMaximumSupportedRefinementPrecision)) {}

  bool addSegment(double first_x, double first_y, double second_x, double second_y) {
    if (!valid_ || !std::isfinite(first_x) || !std::isfinite(first_y) ||
        !std::isfinite(second_x) || !std::isfinite(second_y)) {
      valid_ = false;
      return false;
    }

    const ExactDyadic dx = subtract(exactDyadic(second_x), exactDyadic(first_x));
    const ExactDyadic dy = subtract(exactDyadic(second_y), exactDyadic(first_y));
    const ExactDyadic squared_length = add(multiply(dx, dx), multiply(dy, dy));
    if (squared_length.significand == 0)
      return true;

    double candidate = std::hypot(second_x - first_x, second_y - first_y);
    if (!std::isfinite(candidate)) {
      valid_ = false;
      return false;
    }
    if (candidate == 0.0)
      candidate = std::numeric_limits<double>::denorm_min();

    // A high-quality platform hypot is already adjacent to the exact value.
    // Normalize both upward and downward so the returned value is the unique
    // smallest binary64 number whose exact square covers the exact dyadic
    // squared length.  If libc starts farther away, a fixed-size bit-space
    // binary search preserves both bounded work and cross-libc determinism.
    constexpr unsigned int kNearbyCandidateAttempts = 8;
    if (candidateSquareCovers(candidate, squared_length)) {
      for (unsigned int attempt = 0; attempt < kNearbyCandidateAttempts; ++attempt) {
        const double previous =
          std::nextafter(candidate, -std::numeric_limits<double>::infinity());
        if (previous < 0.0 || !candidateSquareCovers(previous, squared_length))
          break;
        candidate = previous;
      }
      const double previous =
        std::nextafter(candidate, -std::numeric_limits<double>::infinity());
      if (previous >= 0.0 && candidateSquareCovers(previous, squared_length) &&
          !smallestSegmentUpperBound(squared_length, candidate)) {
        valid_ = false;
        return false;
      }
    } else {
      bool repaired = false;
      for (unsigned int attempt = 0; attempt < kNearbyCandidateAttempts; ++attempt) {
        candidate = std::nextafter(candidate, std::numeric_limits<double>::infinity());
        if (candidateSquareCovers(candidate, squared_length)) {
          repaired = true;
          break;
        }
      }
      if (!repaired && !smallestSegmentUpperBound(squared_length, candidate)) {
        valid_ = false;
        return false;
      }
    }

    const ExactDyadic segment_upper = exactDyadic(candidate);
    accumulated_upper_bound_ = add(accumulated_upper_bound_, segment_upper);
    const double previous =
      std::nextafter(candidate, -std::numeric_limits<double>::infinity());
    if (previous < 0.0 || candidateSquareCovers(previous, squared_length)) {
      valid_ = false;
      return false;
    }
    ExactDyadic exact_root;
    if (exactDyadicSquareRoot(squared_length, exact_root)) {
      accumulated_exact_length_ = add(accumulated_exact_length_, exact_root);
      accumulated_lower_bound_ = add(
        accumulated_lower_bound_,
        compare(segment_upper, exact_root) == 0 ? segment_upper : exactDyadic(previous));
    } else {
      accumulated_lower_bound_ = add(accumulated_lower_bound_, exactDyadic(previous));
      irrational_squared_lengths_.emplace_back(squared_length);
    }
    return true;
  }

  // Return a binary64 value which is proved no greater than the mathematical
  // polyline length. It need not be the greatest such value; keeping this
  // operation cheap is important because the planner evaluates every
  // frontier candidate with it.
  bool lowerBound(double& value) const {
    value = std::numeric_limits<double>::quiet_NaN();
    if (!valid_)
      return false;
    if (accumulated_lower_bound_.significand == 0) {
      value = 0.0;
      return true;
    }
    if (accumulated_lower_bound_.significand < 0)
      return false;
    return floorPositiveDyadic(accumulated_lower_bound_, value);
  }

  bool upperBound(double& value,
                  const std::function<bool()>& stop_requested = {}) const {
    value = std::numeric_limits<double>::quiet_NaN();
    if (!valid_)
      return false;
    if (accumulated_upper_bound_.significand == 0) {
      value = 0.0;
      return true;
    }
    if (accumulated_upper_bound_.significand < 0)
      return false;

    // When every radical has an exact dyadic square root, direct conversion
    // is both sufficient and exact. Otherwise use the old per-segment-ceiling
    // sum only as a proved finite search bracket, then binary-search for the
    // first binary64 candidate covering the complete radical sum.
    if (irrational_squared_lengths_.empty())
      return ceilingPositiveDyadic(accumulated_exact_length_, value);

    double loose_lower = 0.0;
    double loose_upper = 0.0;
    if (!floorPositiveDyadic(accumulated_lower_bound_, loose_lower))
      return false;
    if (!ceilingPositiveDyadic(accumulated_upper_bound_, loose_upper)) {
      loose_upper = std::numeric_limits<double>::max();
      Comparison maximum_comparison = Comparison::equal;
      if (!compareTo(loose_upper, maximum_comparison, stop_requested) ||
          maximum_comparison == Comparison::greater) {
        return false;
      }
    }
    std::uint64_t lower_bits = binary64Bits(loose_lower);
    std::uint64_t upper_bits = binary64Bits(loose_upper);
    if (lower_bits > upper_bits)
      return false;
    while (lower_bits < upper_bits) {
      if (stop_requested && stop_requested())
        return false;
      const std::uint64_t middle = lower_bits + (upper_bits - lower_bits) / 2;
      Comparison comparison = Comparison::equal;
      if (!compareTo(binary64FromBits(middle), comparison, stop_requested))
        return false;
      if (comparison == Comparison::greater)
        lower_bits = middle + 1;
      else
        upper_bits = middle;
    }
    value = binary64FromBits(lower_bits);
    Comparison final_comparison = Comparison::equal;
    return std::isfinite(value) && compareTo(value, final_comparison, stop_requested) &&
           final_comparison != Comparison::greater;
  }

  // Compare the exact mathematical polyline length with a nonnegative finite
  // binary64 bound. For irrational sums, fixed-point square-root intervals
  // are doubled until disjoint from the dyadic bound. A positive sum of
  // rational square roots can equal a rational number only when all of its
  // square roots are rational, so equality is handled entirely by the exact
  // fast path and cannot make the refinement loop stall.
  bool compareTo(double bound,
                 Comparison& comparison,
                 const std::function<bool()>& stop_requested = {}) const {
    comparison = Comparison::equal;
    if (!valid_ || !std::isfinite(bound) || bound < 0.0)
      return false;
    const ExactDyadic exact_bound = exactDyadic(bound);
    if (irrational_squared_lengths_.empty()) {
      const int order = compare(accumulated_exact_length_, exact_bound);
      comparison = order < 0 ? Comparison::less
                             : (order > 0 ? Comparison::greater : Comparison::equal);
      return true;
    }

    // Binary64 inputs have bounded exponents. 16384 fractional bits is far
    // beyond what is needed to separate adjacent binary64 values even after
    // the worst endpoint scaling; retaining a hard ceiling makes malformed
    // adversarial requests fail closed instead of consuming unbounded work.
    constexpr unsigned int kInitialPrecision = 64;
    if (maximum_refinement_precision_ == 0)
      return false;
    unsigned int precision = std::min(kInitialPrecision, maximum_refinement_precision_);
    while (precision <= maximum_refinement_precision_) {
      if (stop_requested && stop_requested())
        return false;
      ExactDyadic interval_lower;
      ExactDyadic interval_upper;
      if (!fixedPointInterval(
            precision, interval_lower, interval_upper, stop_requested))
        return false;
      if (compare(interval_upper, exact_bound) <= 0) {
        comparison = Comparison::less;
        return true;
      }
      if (compare(interval_lower, exact_bound) > 0) {
        comparison = Comparison::greater;
        return true;
      }
      if (precision == maximum_refinement_precision_)
        break;
      const unsigned int remaining = maximum_refinement_precision_ - precision;
      precision += std::min(precision, remaining);
    }
    return false;
  }

private:
  // Binary64 endpoint exponents make far less precision sufficient in
  // practice. Enforce the documented hard work ceiling even for arbitrary
  // direct-Core constructor inputs: besides bounding cpp_int allocation, it
  // keeps every fixed-point exponent calculation within signed int.
  static constexpr unsigned int kMaximumSupportedRefinementPrecision = 16384;
  static constexpr int kMaximumBinary64SquaredDyadicExponent =
    2 * (std::numeric_limits<double>::max_exponent -
         std::numeric_limits<double>::digits);
  static_assert(
    kMaximumSupportedRefinementPrecision <=
      static_cast<unsigned int>(
        (std::numeric_limits<int>::max() - kMaximumBinary64SquaredDyadicExponent) / 2),
    "refinement precision must keep squared dyadic exponents in signed-int range");

  struct ExactDyadic {
    boost::multiprecision::cpp_int significand{0};
    int exponent = 0;
  };

  static ExactDyadic exactDyadic(double value) {
    if (value == 0.0)
      return {};
    int exponent = 0;
    const double fraction = std::frexp(value, &exponent);
    constexpr int digits = std::numeric_limits<double>::digits;
    const double scaled = std::ldexp(fraction, digits);
    return ExactDyadic{
      boost::multiprecision::cpp_int(static_cast<std::int64_t>(scaled)), exponent - digits};
  }

  static ExactDyadic add(const ExactDyadic& lhs, const ExactDyadic& rhs) {
    if (lhs.significand == 0)
      return rhs;
    if (rhs.significand == 0)
      return lhs;
    const int exponent = std::min(lhs.exponent, rhs.exponent);
    return ExactDyadic{
      (lhs.significand << (lhs.exponent - exponent)) +
        (rhs.significand << (rhs.exponent - exponent)),
      exponent};
  }

  static ExactDyadic subtract(const ExactDyadic& lhs, const ExactDyadic& rhs) {
    if (lhs.significand == 0)
      return ExactDyadic{-rhs.significand, rhs.exponent};
    if (rhs.significand == 0)
      return lhs;
    const int exponent = std::min(lhs.exponent, rhs.exponent);
    return ExactDyadic{
      (lhs.significand << (lhs.exponent - exponent)) -
        (rhs.significand << (rhs.exponent - exponent)),
      exponent};
  }

  static ExactDyadic multiply(const ExactDyadic& lhs, const ExactDyadic& rhs) {
    if (lhs.significand == 0 || rhs.significand == 0)
      return {};
    return ExactDyadic{lhs.significand * rhs.significand, lhs.exponent + rhs.exponent};
  }

  static bool exactDyadicSquareRoot(const ExactDyadic& squared_length,
                                    ExactDyadic& root) {
    root = {};
    if (squared_length.significand < 0)
      return false;
    if (squared_length.significand == 0)
      return true;
    boost::multiprecision::cpp_int adjusted_significand = squared_length.significand;
    int adjusted_exponent = squared_length.exponent;
    if (adjusted_exponent % 2 != 0) {
      adjusted_significand <<= 1;
      --adjusted_exponent;
    }
    const boost::multiprecision::cpp_int candidate =
      boost::multiprecision::sqrt(adjusted_significand);
    if (candidate * candidate != adjusted_significand)
      return false;
    root = ExactDyadic{candidate, adjusted_exponent / 2};
    return true;
  }

  static int compare(const ExactDyadic& lhs, const ExactDyadic& rhs) {
    if (lhs.significand == 0 && rhs.significand == 0)
      return 0;
    const int exponent = std::min(lhs.exponent, rhs.exponent);
    const boost::multiprecision::cpp_int lhs_aligned =
      lhs.significand << (lhs.exponent - exponent);
    const boost::multiprecision::cpp_int rhs_aligned =
      rhs.significand << (rhs.exponent - exponent);
    if (lhs_aligned < rhs_aligned)
      return -1;
    if (lhs_aligned > rhs_aligned)
      return 1;
    return 0;
  }

  static bool candidateSquareCovers(double candidate, const ExactDyadic& squared_length) {
    if (!std::isfinite(candidate) || candidate < 0.0)
      return false;
    const ExactDyadic exact_candidate = exactDyadic(candidate);
    return compare(multiply(exact_candidate, exact_candidate), squared_length) >= 0;
  }

  static double binary64FromBits(std::uint64_t bits) {
    double value = 0.0;
    static_assert(sizeof(value) == sizeof(bits), "binary64 storage must be 64 bits");
    std::memcpy(&value, &bits, sizeof(value));
    return value;
  }

  static std::uint64_t binary64Bits(double value) {
    std::uint64_t bits = 0;
    static_assert(sizeof(value) == sizeof(bits), "binary64 storage must be 64 bits");
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
  }

  static bool smallestSegmentUpperBound(const ExactDyadic& squared_length, double& result) {
    constexpr std::uint64_t kMaximumFinitePositiveBits = UINT64_C(0x7fefffffffffffff);
    const double maximum = binary64FromBits(kMaximumFinitePositiveBits);
    if (!candidateSquareCovers(maximum, squared_length))
      return false;

    // For nonnegative binary64 values, the raw bit pattern has the same total
    // order as the represented value.  Keep zero as the proven-insufficient
    // lower endpoint because this helper is called only for nonzero lengths.
    std::uint64_t lower = 0;
    std::uint64_t upper = kMaximumFinitePositiveBits;
    while (lower + 1 < upper) {
      const std::uint64_t middle = lower + (upper - lower) / 2;
      if (candidateSquareCovers(binary64FromBits(middle), squared_length))
        upper = middle;
      else
        lower = middle;
    }
    result = binary64FromBits(upper);
    return true;
  }

  static boost::multiprecision::cpp_int ceilDivideByPowerOfTwo(
    const boost::multiprecision::cpp_int& value, unsigned int shift) {
    if (shift == 0)
      return value;
    boost::multiprecision::cpp_int quotient = value >> shift;
    if ((quotient << shift) != value)
      ++quotient;
    return quotient;
  }

  static boost::multiprecision::cpp_int floorScaledDyadic(const ExactDyadic& value,
                                                           unsigned int precision) {
    const int scaled_exponent = value.exponent + static_cast<int>(precision);
    if (scaled_exponent >= 0)
      return value.significand << scaled_exponent;
    return value.significand >> static_cast<unsigned int>(-scaled_exponent);
  }

  static boost::multiprecision::cpp_int ceilScaledDyadic(const ExactDyadic& value,
                                                          unsigned int precision) {
    const int scaled_exponent = value.exponent + static_cast<int>(precision);
    if (scaled_exponent >= 0)
      return value.significand << scaled_exponent;
    return ceilDivideByPowerOfTwo(
      value.significand, static_cast<unsigned int>(-scaled_exponent));
  }

  static boost::multiprecision::cpp_int scaledSquareRootFloor(
    const ExactDyadic& squared_length, unsigned int precision) {
    const int scaled_exponent =
      squared_length.exponent + 2 * static_cast<int>(precision);
    boost::multiprecision::cpp_int scaled_radicand = squared_length.significand;
    if (scaled_exponent >= 0)
      scaled_radicand <<= scaled_exponent;
    else
      scaled_radicand >>= static_cast<unsigned int>(-scaled_exponent);
    return boost::multiprecision::sqrt(scaled_radicand);
  }

  bool fixedPointInterval(unsigned int precision,
                          ExactDyadic& interval_lower,
                          ExactDyadic& interval_upper,
                          const std::function<bool()>& stop_requested) const {
    if (!valid_ || precision > static_cast<unsigned int>(std::numeric_limits<int>::max()))
      return false;
    boost::multiprecision::cpp_int lower =
      floorScaledDyadic(accumulated_exact_length_, precision);
    boost::multiprecision::cpp_int upper =
      ceilScaledDyadic(accumulated_exact_length_, precision);
    for (const auto& squared_length : irrational_squared_lengths_) {
      if (stop_requested && stop_requested())
        return false;
      const boost::multiprecision::cpp_int segment_lower =
        scaledSquareRootFloor(squared_length, precision);
      lower += segment_lower;
      // Irrational roots cannot lie on a dyadic fixed-point grid, hence the
      // next integer is a strict upper bound at every precision.
      upper += segment_lower + 1;
    }
    interval_lower = ExactDyadic{std::move(lower), -static_cast<int>(precision)};
    interval_upper = ExactDyadic{std::move(upper), -static_cast<int>(precision)};
    return true;
  }

  static bool ceilingPositiveDyadic(const ExactDyadic& exact_value, double& result) {
    if (exact_value.significand <= 0)
      return false;

    const int most_significant_bit =
      static_cast<int>(boost::multiprecision::msb(exact_value.significand));
    const int binary_exponent = exact_value.exponent + most_significant_bit;
    if (binary_exponent > std::numeric_limits<double>::max_exponent - 1)
      return false;

    constexpr int kMinimumNormalExponent = std::numeric_limits<double>::min_exponent - 1;
    constexpr int kSubnormalQuantumExponent =
      std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits;
    const int quantum_exponent =
      binary_exponent >= kMinimumNormalExponent
        ? binary_exponent - (std::numeric_limits<double>::digits - 1)
        : kSubnormalQuantumExponent;

    boost::multiprecision::cpp_int rounded_significand;
    if (exact_value.exponent >= quantum_exponent) {
      rounded_significand =
        exact_value.significand << (exact_value.exponent - quantum_exponent);
    } else {
      rounded_significand = ceilDivideByPowerOfTwo(
        exact_value.significand,
        static_cast<unsigned int>(quantum_exponent - exact_value.exponent));
    }

    const boost::multiprecision::cpp_int maximum_rounded_significand =
      boost::multiprecision::cpp_int(1) << std::numeric_limits<double>::digits;
    if (rounded_significand <= 0 || rounded_significand > maximum_rounded_significand)
      return false;

    const double significand = rounded_significand.convert_to<double>();
    const double candidate = std::ldexp(significand, quantum_exponent);
    if (!std::isfinite(candidate))
      return false;
    result = candidate;
    return compare(exactDyadic(result), exact_value) >= 0;
  }

  static bool floorPositiveDyadic(const ExactDyadic& exact_value, double& result) {
    if (exact_value.significand < 0)
      return false;
    if (exact_value.significand == 0) {
      result = 0.0;
      return true;
    }

    const int most_significant_bit =
      static_cast<int>(boost::multiprecision::msb(exact_value.significand));
    const int binary_exponent = exact_value.exponent + most_significant_bit;
    if (binary_exponent > std::numeric_limits<double>::max_exponent - 1) {
      result = std::numeric_limits<double>::max();
      return compare(exactDyadic(result), exact_value) <= 0;
    }

    constexpr int kMinimumNormalExponent = std::numeric_limits<double>::min_exponent - 1;
    constexpr int kSubnormalQuantumExponent =
      std::numeric_limits<double>::min_exponent - std::numeric_limits<double>::digits;
    const int quantum_exponent =
      binary_exponent >= kMinimumNormalExponent
        ? binary_exponent - (std::numeric_limits<double>::digits - 1)
        : kSubnormalQuantumExponent;

    boost::multiprecision::cpp_int rounded_significand;
    if (exact_value.exponent >= quantum_exponent) {
      rounded_significand =
        exact_value.significand << (exact_value.exponent - quantum_exponent);
    } else {
      rounded_significand =
        exact_value.significand >> (quantum_exponent - exact_value.exponent);
    }
    if (rounded_significand == 0) {
      result = 0.0;
      return true;
    }

    const boost::multiprecision::cpp_int maximum_rounded_significand =
      boost::multiprecision::cpp_int(1) << std::numeric_limits<double>::digits;
    if (rounded_significand < 0 || rounded_significand > maximum_rounded_significand)
      return false;
    const double significand = rounded_significand.convert_to<double>();
    const double candidate = std::ldexp(significand, quantum_exponent);
    if (!std::isfinite(candidate))
      return false;
    result = candidate;
    return compare(exactDyadic(result), exact_value) <= 0;
  }

  bool valid_{true};
  ExactDyadic accumulated_lower_bound_;
  ExactDyadic accumulated_upper_bound_;
  ExactDyadic accumulated_exact_length_;
  std::vector<ExactDyadic> irrational_squared_lengths_;
  unsigned int maximum_refinement_precision_;
};

}  // namespace raystar
