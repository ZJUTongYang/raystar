#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_2.h>

#include <optional>

namespace raystar
{
namespace exact_geometry
{

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = Kernel::Point_2;
using FT = Kernel::FT;

// Intersect a finite blocker segment with the forward ray from `source`
// through `ray_anchor`. The old visibility code solved the two infinite line
// equations; keeping the parameter domains here makes the same construction
// exact while rejecting an invalid hit behind the source or beyond the
// obstacle edge. Collinear overlap is ambiguous for a single-point endpoint
// and is therefore reported as an ordinary failure.
inline std::optional<Point> intersectSegmentWithRay(
  const Point& segment_start, const Point& segment_end,
  const Point& source, const Point& ray_anchor)
{
  if (segment_start == segment_end || source == ray_anchor)
    return std::nullopt;

  const FT a1 = segment_end.y() - segment_start.y();
  const FT b1 = -(segment_end.x() - segment_start.x());
  const FT c1 = -a1 * segment_start.x() - b1 * segment_start.y();
  const FT a2 = ray_anchor.y() - source.y();
  const FT b2 = -(ray_anchor.x() - source.x());
  const FT c2 = -a2 * source.x() - b2 * source.y();
  const FT determinant = a1 * b2 - a2 * b1;
  if (determinant == FT(0))
    return std::nullopt;

  const Point intersection(
    (b1 * c2 - b2 * c1) / determinant,
    (c1 * a2 - c2 * a1) / determinant);
  if (!Kernel::Segment_2(segment_start, segment_end).has_on(intersection) ||
      !Kernel::Ray_2(source, ray_anchor).has_on(intersection))
  {
    return std::nullopt;
  }
  return intersection;
}

enum class PortalRayPosition
{
  before_lower,
  equal_lower,
  inside,
  equal_upper,
  after_upper,
  invalid
};

inline bool isSameDirectedRay(const Point& source, const Point& first,
  const Point& second)
{
  if (first == source || second == source)
    return false;
  return CGAL::orientation(source, first, second) == CGAL::COLLINEAR &&
    CGAL::scalar_product(first - source, second - source) > FT(0);
}

// Compare directions by their counter-clockwise offset from reference in
// [0, 2*pi). Distance from source is intentionally ignored: points on the
// same directed ray compare equal so stable_sort can preserve visibility
// discontinuities (far point followed by near point, or vice versa).
inline CGAL::Comparison_result compareRaysCounterClockwiseFrom(
  const Point& source, const Point& reference,
  const Point& first, const Point& second)
{
  if (first == second)
    return CGAL::EQUAL;
  if (first == source)
    return second == source ? CGAL::EQUAL : CGAL::SMALLER;
  if (second == source)
    return CGAL::LARGER;
  if (reference == source)
    return CGAL::EQUAL;

  const auto half = [&](const Point& point) {
    const CGAL::Orientation side = CGAL::orientation(source, reference, point);
    if (side == CGAL::LEFT_TURN)
      return 0;
    if (side == CGAL::RIGHT_TURN)
      return 1;
    return CGAL::scalar_product(reference - source, point - source) >= FT(0) ? 0 : 1;
  };

  const int first_half = half(first);
  const int second_half = half(second);
  if (first_half != second_half)
    return first_half < second_half ? CGAL::SMALLER : CGAL::LARGER;

  const CGAL::Orientation turn = CGAL::orientation(source, first, second);
  if (turn == CGAL::LEFT_TURN)
    return CGAL::SMALLER;
  if (turn == CGAL::RIGHT_TURN)
    return CGAL::LARGER;

  // Within one half-plane, collinear non-zero vectors have the same direction.
  return CGAL::EQUAL;
}

// Exact ordering equivalent to sorting std::atan2(y, x) in ascending order:
// (-pi, 0) first, then [0, pi]. In particular the negative x-axis (atan2=pi)
// remains at the end, matching the former implementation.
inline CGAL::Comparison_result compareRaysLikeAtan2(
  const Point& source, const Point& first, const Point& second)
{
  if (first == second)
    return CGAL::EQUAL;
  if (first == source)
    return second == source ? CGAL::EQUAL : CGAL::SMALLER;
  if (second == source)
    return CGAL::LARGER;

  const auto lower_half = [&](const Point& point) {
    return (point - source).y() < FT(0);
  };
  const bool first_lower = lower_half(first);
  const bool second_lower = lower_half(second);
  if (first_lower != second_lower)
    return first_lower ? CGAL::SMALLER : CGAL::LARGER;

  const CGAL::Orientation turn = CGAL::orientation(source, first, second);
  if (turn == CGAL::LEFT_TURN)
    return CGAL::SMALLER;
  if (turn == CGAL::RIGHT_TURN)
    return CGAL::LARGER;

  const FT dot = CGAL::scalar_product(first - source, second - source);
  if (dot > FT(0))
    return CGAL::EQUAL;

  // The only opposite directions in the same atan2 half are the two x-axis
  // directions. Positive x (angle 0) precedes negative x (angle pi).
  return (first - source).x() > FT(0) ? CGAL::SMALLER : CGAL::LARGER;
}

inline bool isClosedCounterClockwiseSweepMember(
  const Point& source, const Point& start, const Point& end,
  const Point& candidate)
{
  if (start == source || end == source || candidate == source)
    return false;
  return compareRaysCounterClockwiseFrom(source, start, candidate, end) != CGAL::LARGER;
}

// Reproduce the former Bungiu CASE dispatch exactly, but with exact predicates.
// The algorithm's portal invariant is a counter-clockwise span in [0, pi].
// A clockwise upper ray would have been represented as a span greater than pi
// while theta was restricted to (lower-pi, lower+pi], so it is rejected.
inline PortalRayPosition classifyPortalRay(
  const Point& source, const Point& lower, const Point& upper,
  const Point& candidate)
{
  if (lower == source || upper == source || candidate == source)
    return PortalRayPosition::invalid;

  const CGAL::Orientation upper_turn = CGAL::orientation(source, lower, upper);
  if (upper_turn == CGAL::RIGHT_TURN)
    return PortalRayPosition::invalid;

  const bool zero_span = isSameDirectedRay(source, lower, upper);
  const bool pi_span = upper_turn == CGAL::COLLINEAR && !zero_span;

  if (isSameDirectedRay(source, lower, candidate))
    return PortalRayPosition::equal_lower;

  const CGAL::Orientation candidate_turn =
    CGAL::orientation(source, lower, candidate);
  if (candidate_turn == CGAL::RIGHT_TURN)
    return PortalRayPosition::before_lower;

  if (!zero_span && isSameDirectedRay(source, upper, candidate))
    return PortalRayPosition::equal_upper;

  if (zero_span)
    return PortalRayPosition::after_upper;
  if (pi_span)
    return PortalRayPosition::inside;

  return CGAL::orientation(source, candidate, upper) == CGAL::LEFT_TURN ?
    PortalRayPosition::inside : PortalRayPosition::after_upper;
}

inline bool isPortalSpanAtMostPi(
  const Point& source, const Point& lower, const Point& upper)
{
  if (lower == source || upper == source)
    return false;
  return CGAL::orientation(source, lower, upper) != CGAL::RIGHT_TURN;
}

inline bool isRemovableCollinearMiddle(
  const Point& previous, const Point& current, const Point& next)
{
  return CGAL::orientation(previous, current, next) == CGAL::COLLINEAR &&
    CGAL::collinear_are_ordered_along_line(previous, current, next);
}

}  // namespace exact_geometry
}  // namespace raystar
