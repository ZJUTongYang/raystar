# Result for one entry of a multi-goal bounded request. goal_index and goal
# preserve the input association even when goals are duplicated.
uint32                                goal_index
geometry_msgs/PoseStamped             goal
float64                               requested_max_path_length
# Payload presence only: true iff path_results is non-empty. A false value can
# still accompany a completely satisfied bounded enumeration; use result_info
# for the per-goal search/output certificate.
bool                                  success
raystar_interfaces/PlanningResultInfo result_info
string                                message
raystar_interfaces/PathResult[]       path_results
