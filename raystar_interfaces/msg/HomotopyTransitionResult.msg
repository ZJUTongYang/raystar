uint8 STATUS_SUCCESS=0
uint8 STATUS_INVALID_REFERENCE=1
uint8 STATUS_NO_CORRIDOR=2
uint8 STATUS_STOPPED=3
uint8 STATUS_FAILURE=4

raystar_interfaces/ConfigurationTransitionPair pair
uint8 status
nav_msgs/Path path
# For STATUS_SUCCESS, the least binary64 value that conservatively upper-bounds
# the mathematical Euclidean length of the serialized world-coordinate `path`.
# A successful one-pose identity transition has length zero. Consumers must
# not use this field as an edge cost for any non-success status.
float64 path_length
bool collision_free
bool homotopy_preserved
bool locally_shortest

# Stable CDT face IDs in lifted-sleeve occurrence order. IDs may repeat when
# a transition winds through the same geometric triangle more than once.
uint32[] triangle_occurrences
string message
