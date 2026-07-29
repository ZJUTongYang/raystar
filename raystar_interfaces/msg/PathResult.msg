# One planned path and its metric cost. The path and cost are deliberately one
# record so they cannot become misaligned during filtering or serialization.
nav_msgs/Path path
float64       cost  # Path length in metres
