#pragma once

#include <optional>
#include <string>

#include <raystar/cooperative_stop.h>
#include <raystar/polymap.h>

namespace raystar
{
namespace detail
{

// A visibility boundary has two different closure semantics.  A free-space
// source produces a closed cycle, while an obstacle-vertex/scoped source
// produces an open angular sector whose conceptual closing edges pass through
// the source.
enum class VisibilityBoundaryMode
{
  closed_cycle,
  open_sector
};

struct VisibilityGeometryContext
{
  exact_geometry::Point source;
  VisibilityBoundaryMode mode;
  std::optional<exact_geometry::Point> start_anchor;
  std::optional<exact_geometry::Point> end_anchor;
};

// Validate only source-dependent geometry.  Boundary support/topology metadata
// remains the responsibility of Polymap::validateVisibilityRegion().  The
// accepted boundary is intentionally weak: repeated endpoints, endpoint
// touches and collinear retracing are legal visibility discontinuities.  A
// strict interior crossing, however, is never legal.  This is the generated
// visibility invariant used by Raystar, not a general-purpose formal
// weak-simplicity recognizer.
bool validateVisibilityGeometry(
  const VisibilityRegion& visibility_region,
  const VisibilityGeometryContext& context,
  std::string* error = nullptr);

// Cooperative variant.  A stopped operation never returns partially
// validated geometry and leaves error empty; ordinary validation failures keep
// the same diagnostic messages as the bool wrapper above.
OperationStatus validateVisibilityGeometry(
  const VisibilityRegion& visibility_region,
  const VisibilityGeometryContext& context,
  const StopToken& stop_token,
  std::string* error = nullptr);

}  // namespace detail
}  // namespace raystar
