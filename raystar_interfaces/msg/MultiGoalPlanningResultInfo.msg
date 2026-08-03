# Aggregate status for one shared-tree multi-goal request. status and
# limits_reached use the constants defined by PlanningResultInfo. Per-goal
# completion and output certificates live in GoalPathResult.result_info.
unique_identifier_msgs/UUID map_id
unique_identifier_msgs/UUID environment_id
uint8  status
uint16 limits_reached
bool   request_satisfied     # Every goal search and bounded output is complete, paths optional
bool   search_complete
bool   output_complete       # Every per-goal authoritative metric-eligible path was returned
bool   debug_requested
bool   debug_output_complete
uint32 requested_goal_count
uint32 returned_goal_count
uint32 completed_goal_count
uint32 goals_with_paths
uint32 found_path_count      # Sum of per-goal metric-eligible paths before resource filtering
uint32 returned_path_count
uint64 expanded_nodes
float64 map_time_ms
float64 plan_time_ms
