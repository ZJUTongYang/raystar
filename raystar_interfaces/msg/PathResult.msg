# One planned path and its metric cost. The path and cost are deliberately one
# record so they cannot become misaligned during filtering or serialization.
#
# path is the densely sampled consumer path used for execution and display.
# topology_path preserves the unsampled piecewise-linear Core representation
# (continuous endpoints plus obstacle turning vertices). cost is the least
# binary64 value that conservatively upper-bounds the mathematical polyline
# length of both serialized representations; interpolation/world-coordinate
# rounding can make their exact lengths differ microscopically. Algorithms
# that consume a path as a homotopy reference, including UPS, must use
# topology_path so interpolation cannot perturb a boundary-following segment.
# Both paths are oriented from the request start to its goal. Each returned
# path-results sequence is sorted by nondecreasing final public cost.
nav_msgs/Path path
nav_msgs/Path topology_path
float64       cost  # Conservative serialized-path length ceiling in metres
