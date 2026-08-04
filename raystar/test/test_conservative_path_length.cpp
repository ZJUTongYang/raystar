#include "conservative_path_length.h"

#include <gtest/gtest.h>

#include <cstddef>
#include <cfenv>
#include <cmath>
#include <limits>

namespace raystar {
namespace {

TEST(ConservativeBinary64PathLengthTest, CatchesHypotRoundingBelowInclusiveBound) {
  // The serialized endpoint is exact binary64, so the squared segment length
  // is exactly 1 + 2^-54.  Its length is strictly greater than 1, while the
  // nearest binary64 hypot rounds down to 1 and used to pass a bound of 1.
  const double dy = std::ldexp(1.0, -27);
  ASSERT_DOUBLE_EQ(std::hypot(1.0, dy), 1.0);

  ConservativeBinary64PathLength certificate;
  ASSERT_TRUE(certificate.addSegment(0.0, 0.0, 1.0, dy));
  double upper_bound = std::numeric_limits<double>::quiet_NaN();
  ASSERT_TRUE(certificate.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound, std::nextafter(1.0, std::numeric_limits<double>::infinity()));
  EXPECT_GT(upper_bound, 1.0);

  // sqrt(2) rounds upward on the test platform.  The same exact-square
  // normalization must retain that smallest covering value, independent of
  // whether a different libc starts one ULP higher or lower.
  ConservativeBinary64PathLength square_root_two;
  ASSERT_TRUE(square_root_two.addSegment(0.0, 0.0, 1.0, 1.0));
  ASSERT_TRUE(square_root_two.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound, 0x1.6a09e667f3bcdp+0);
}

TEST(ConservativeBinary64PathLengthTest, RetainsExactAxisAndPythagoreanLengths) {
  ConservativeBinary64PathLength certificate;
  ASSERT_TRUE(certificate.addSegment(0.0, 0.0, 3.0, 4.0));
  ASSERT_TRUE(certificate.addSegment(3.0, 4.0, 6.0, 8.0));
  double upper_bound = std::numeric_limits<double>::quiet_NaN();
  ASSERT_TRUE(certificate.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound, 10.0);

  ConservativeBinary64PathLength subnormal;
  const double quantum = std::numeric_limits<double>::denorm_min();
  ASSERT_TRUE(subnormal.addSegment(0.0, 0.0, quantum, 0.0));
  ASSERT_TRUE(subnormal.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound, quantum);
}

TEST(ConservativeBinary64PathLengthTest, RefinesTheCompleteRadicalSumAtOneUlpBoundary) {
  // The exact length is (sqrt(2026) + sqrt(5)) / 2. Both individually
  // rounded segment ceilings sum to the next binary64 value, so accumulating
  // per-segment ceilings would incorrectly reject this inclusive bound.
  constexpr double inclusive_bound = 0x1.79fa384f9da53p+4;
  ConservativeBinary64PathLength certificate;
  ASSERT_TRUE(certificate.addSegment(19.5, 52.5, 20.0, 30.0));
  ASSERT_TRUE(certificate.addSegment(20.0, 30.0, 20.5, 29.0));

  double lower_bound = std::numeric_limits<double>::quiet_NaN();
  double upper_bound = std::numeric_limits<double>::quiet_NaN();
  ASSERT_TRUE(certificate.lowerBound(lower_bound));
  ASSERT_TRUE(certificate.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(lower_bound, 0x1.79fa384f9da52p+4);
  EXPECT_DOUBLE_EQ(upper_bound, inclusive_bound);

  ConservativeBinary64PathLength::Comparison comparison =
    ConservativeBinary64PathLength::Comparison::equal;
  ASSERT_TRUE(certificate.compareTo(inclusive_bound, comparison));
  EXPECT_EQ(comparison, ConservativeBinary64PathLength::Comparison::less);
  ASSERT_TRUE(certificate.compareTo(std::nextafter(inclusive_bound, 0.0), comparison));
  EXPECT_EQ(comparison, ConservativeBinary64PathLength::Comparison::greater);
}

TEST(ConservativeBinary64PathLengthTest, TightBoundsIgnoreProcessRoundingMode) {
  const int original_rounding = std::fegetround();
  ASSERT_NE(original_rounding, -1);
  struct RoundingModeGuard {
    int mode;
    ~RoundingModeGuard() {
      (void)std::fesetround(mode);
    }
  } rounding_guard{original_rounding};
  (void)rounding_guard;
  constexpr double expected_upper = 0x1.79fa384f9da53p+4;
  for (const int rounding_mode : {FE_DOWNWARD, FE_UPWARD, FE_TOWARDZERO, FE_TONEAREST}) {
    ASSERT_EQ(std::fesetround(rounding_mode), 0);
    ConservativeBinary64PathLength certificate;
    ASSERT_TRUE(certificate.addSegment(19.5, 52.5, 20.0, 30.0));
    ASSERT_TRUE(certificate.addSegment(20.0, 30.0, 20.5, 29.0));
    double lower_bound = 0.0;
    double upper_bound = 0.0;
    ASSERT_TRUE(certificate.lowerBound(lower_bound));
    ASSERT_TRUE(certificate.upperBound(upper_bound));
    EXPECT_DOUBLE_EQ(lower_bound, 0x1.79fa384f9da52p+4);
    EXPECT_DOUBLE_EQ(upper_bound, expected_upper);
  }
}

TEST(ConservativeBinary64PathLengthTest, ExactComparisonRejectsLowerBoundFalsePositive) {
  const double dy = std::ldexp(1.0, -27);
  ConservativeBinary64PathLength certificate;
  ASSERT_TRUE(certificate.addSegment(0.0, 0.0, 1.0, dy));

  double lower_bound = std::numeric_limits<double>::quiet_NaN();
  ASSERT_TRUE(certificate.lowerBound(lower_bound));
  EXPECT_DOUBLE_EQ(lower_bound, 1.0);
  ConservativeBinary64PathLength::Comparison comparison =
    ConservativeBinary64PathLength::Comparison::equal;
  ASSERT_TRUE(certificate.compareTo(1.0, comparison));
  EXPECT_EQ(comparison, ConservativeBinary64PathLength::Comparison::greater);
}

TEST(ConservativeBinary64PathLengthTest, UnresolvedAndStoppedRefinementFailClosed) {
  // sqrt(1 + 2^-200) differs from 1 below a 64-bit fixed-point quantum.
  // An intentionally restricted test seam must report unresolved rather than
  // guessing equality or allowing the path through an inclusive bound.
  ConservativeBinary64PathLength restricted(64);
  ASSERT_TRUE(restricted.addSegment(0.0, 0.0, 1.0, std::ldexp(1.0, -100)));
  ConservativeBinary64PathLength::Comparison comparison =
    ConservativeBinary64PathLength::Comparison::equal;
  EXPECT_FALSE(restricted.compareTo(1.0, comparison));

  ConservativeBinary64PathLength stopped;
  ASSERT_TRUE(stopped.addSegment(19.5, 52.5, 20.0, 30.0));
  ASSERT_TRUE(stopped.addSegment(20.0, 30.0, 20.5, 29.0));
  size_t polls = 0;
  const auto request_stop = [&]() { return ++polls >= 2; };
  double upper_bound = 0.0;
  EXPECT_FALSE(stopped.upperBound(upper_bound, request_stop));
  EXPECT_GE(polls, 2u);
}

TEST(ConservativeBinary64PathLengthTest, ArbitraryPrecisionInputIsHardCappedSafely) {
  ConservativeBinary64PathLength certificate(std::numeric_limits<unsigned int>::max());
  ASSERT_TRUE(certificate.addSegment(19.5, 52.5, 20.0, 30.0));
  ASSERT_TRUE(certificate.addSegment(20.0, 30.0, 20.5, 29.0));
  double upper_bound = std::numeric_limits<double>::quiet_NaN();
  ASSERT_TRUE(certificate.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound, 0x1.79fa384f9da53p+4);
  ConservativeBinary64PathLength::Comparison comparison =
    ConservativeBinary64PathLength::Comparison::equal;
  ASSERT_TRUE(certificate.compareTo(std::nextafter(upper_bound, 0.0), comparison));
  EXPECT_EQ(comparison, ConservativeBinary64PathLength::Comparison::greater);
}

TEST(ConservativeBinary64PathLengthTest, RecognizesDyadicRootsWiderThanBinary64) {
  // 1 - denorm_min is an exact dyadic axis length with more than 53
  // significant bits, so it is not itself representable as binary64. It is
  // nevertheless rational and, together with denorm_min, sums exactly to 1.
  const double quantum = std::numeric_limits<double>::denorm_min();
  ConservativeBinary64PathLength certificate;
  ASSERT_TRUE(certificate.addSegment(0.0, 0.0, quantum, 0.0));
  ASSERT_TRUE(certificate.addSegment(quantum, 0.0, 1.0, 0.0));
  double upper_bound = std::numeric_limits<double>::quiet_NaN();
  ASSERT_TRUE(certificate.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound, 1.0);
  ConservativeBinary64PathLength::Comparison comparison =
    ConservativeBinary64PathLength::Comparison::less;
  ASSERT_TRUE(certificate.compareTo(1.0, comparison));
  EXPECT_EQ(comparison, ConservativeBinary64PathLength::Comparison::equal);

  // Exercise the same property with non-axis Pythagorean segments whose
  // exact 54-bit hypotenuses are individually non-representable. Their sum
  // is representable and must take the exact-equality path, not be mistaken
  // for an irrational radical sum.
  constexpr std::int64_t a1 = INT64_C(6075000180000001);
  constexpr std::int64_t b1 = INT64_C(8100000090000000);
  constexpr std::int64_t a2 = INT64_C(5625000300000003);
  constexpr std::int64_t b2 = INT64_C(7500000150000000);
  const double scale = std::ldexp(1.0, -60);
  const double first_x = static_cast<double>(a1) * scale;
  const double first_y = static_cast<double>(b1) * scale;
  const double second_x = static_cast<double>(a1 - a2) * scale;
  const double second_y = static_cast<double>(b1 - b2) * scale;
  ConservativeBinary64PathLength pythagorean;
  ASSERT_TRUE(pythagorean.addSegment(0.0, 0.0, first_x, first_y));
  ASSERT_TRUE(pythagorean.addSegment(first_x, first_y, second_x, second_y));
  ASSERT_TRUE(pythagorean.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound, 0x1.151c96a6ebe01p-6);
  ASSERT_TRUE(pythagorean.compareTo(upper_bound, comparison));
  EXPECT_EQ(comparison, ConservativeBinary64PathLength::Comparison::equal);
}

TEST(ConservativeBinary64PathLengthTest, DirectlyRoundsLargeExactDyadicSumsUpward) {
  ConservativeBinary64PathLength many_small_segments;
  constexpr std::size_t segment_count = 100000;
  const double small_length = std::ldexp(1.0, -40);
  for (std::size_t index = 0; index < segment_count; ++index)
    ASSERT_TRUE(many_small_segments.addSegment(0.0, 0.0, small_length, 0.0));
  double upper_bound = std::numeric_limits<double>::quiet_NaN();
  ASSERT_TRUE(many_small_segments.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound, std::ldexp(static_cast<double>(segment_count), -40));

  // This non-power-of-two value reproduces a 100k-segment case where an
  // 80-bit sequential sum can start more than eight binary64 ULPs below the
  // exact dyadic sum.  Direct conversion has a deterministic exact ceiling.
  ConservativeBinary64PathLength repeated_non_power_of_two;
  constexpr double repeated_length = 0x1.4e8fcf1c6e0d0p+0;
  for (std::size_t index = 0; index < segment_count; ++index) {
    ASSERT_TRUE(repeated_non_power_of_two.addSegment(0.0, 0.0, repeated_length, 0.0));
  }
  ASSERT_TRUE(repeated_non_power_of_two.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound, 0x1.fe802f66c16cap+16);

  // 1024 * 2^-63 is exactly half an ULP above 1.  Directed conversion must
  // choose the next binary64 value instead of ties-to-even rounding down.
  ConservativeBinary64PathLength half_ulp_sum;
  ASSERT_TRUE(half_ulp_sum.addSegment(0.0, 0.0, 1.0, 0.0));
  const double half_ulp_piece = std::ldexp(1.0, -63);
  for (std::size_t index = 0; index < 1024; ++index)
    ASSERT_TRUE(half_ulp_sum.addSegment(0.0, 0.0, half_ulp_piece, 0.0));
  ASSERT_TRUE(half_ulp_sum.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound, std::nextafter(1.0, std::numeric_limits<double>::infinity()));

  // Preserve a tiny exact addend across a 900-bit exponent range.  Any
  // positive addend makes the directed ceiling advance by one large-value ULP.
  ConservativeBinary64PathLength dynamic_range;
  const double large_length = std::ldexp(1.0, 900);
  ASSERT_TRUE(dynamic_range.addSegment(0.0, 0.0, large_length, 0.0));
  ASSERT_TRUE(dynamic_range.addSegment(0.0, 0.0, 1.0, 0.0));
  ASSERT_TRUE(dynamic_range.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound,
                   std::nextafter(large_length, std::numeric_limits<double>::infinity()));
}

TEST(ConservativeBinary64PathLengthTest, FailsClosedWhenNoFiniteBoundExists) {
  ConservativeBinary64PathLength certificate;
  const double maximum = std::numeric_limits<double>::max();
  EXPECT_FALSE(certificate.addSegment(-maximum, 0.0, maximum, 0.0));
  double upper_bound = 0.0;
  EXPECT_FALSE(certificate.upperBound(upper_bound));
  EXPECT_TRUE(std::isnan(upper_bound));

  // Every segment is finite, but their exact dyadic upper-bound sum lies
  // above DBL_MAX and therefore has no finite public binary64 certificate.
  ConservativeBinary64PathLength aggregate_overflow;
  ASSERT_TRUE(aggregate_overflow.addSegment(0.0, 0.0, maximum, 0.0));
  ASSERT_TRUE(
    aggregate_overflow.addSegment(0.0, 0.0, std::numeric_limits<double>::denorm_min(), 0.0));
  upper_bound = 0.0;
  EXPECT_FALSE(aggregate_overflow.upperBound(upper_bound));
  EXPECT_TRUE(std::isnan(upper_bound));
}

TEST(ConservativeBinary64PathLengthTest, TightSumCanFitWhenLooseSegmentCeilingsOverflow) {
  const double maximum = std::numeric_limits<double>::max();
  const double one_third_predecessor = 0x1.5555555555554p+1022;
  ConservativeBinary64PathLength certificate;
  for (int segment = 0; segment < 3; ++segment)
    ASSERT_TRUE(certificate.addSegment(0.0, 0.0, one_third_predecessor, 1.0));

  // Each segment ceiling is the successor of one_third_predecessor, whose
  // three-way dyadic sum exceeds DBL_MAX. The complete radical sum remains
  // below DBL_MAX and has DBL_MAX as its tight binary64 ceiling.
  double upper_bound = std::numeric_limits<double>::quiet_NaN();
  ASSERT_TRUE(certificate.upperBound(upper_bound));
  EXPECT_DOUBLE_EQ(upper_bound, maximum);
  ConservativeBinary64PathLength::Comparison comparison =
    ConservativeBinary64PathLength::Comparison::equal;
  ASSERT_TRUE(certificate.compareTo(maximum, comparison));
  EXPECT_EQ(comparison, ConservativeBinary64PathLength::Comparison::less);

  ConservativeBinary64PathLength true_overflow;
  ASSERT_TRUE(true_overflow.addSegment(0.0, 0.0, maximum, 0.0));
  ASSERT_TRUE(true_overflow.addSegment(0.0, 0.0, std::numeric_limits<double>::denorm_min(), 0.0));
  upper_bound = 0.0;
  EXPECT_FALSE(true_overflow.upperBound(upper_bound));
  EXPECT_TRUE(std::isnan(upper_bound));
}

}  // namespace
}  // namespace raystar
