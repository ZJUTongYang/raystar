"""Contract tests for the wire-facing Raystar ROS interfaces.

The node and RViz plugin are intentionally built against these generated
types.  This test keeps an accidental field rename, reorder, or type change
from silently becoming a wire-compatibility break.  The checks use the
generated Python metadata rather than duplicating the IDL parser, so they also
verify that rosidl generated the expected public schema.
"""

from raystar_interfaces.action import (
    BuildRaystarTransitionGraph,
    PlanRaystarGoalSet,
    PlanRaystarPaths,
)
from raystar_interfaces.msg import (
    ConfigurationTransitionPair,
    DebugNode,
    GoalPathResult,
    HomotopyTransitionResult,
    MapStatus,
    MultiGoalPlanningResultInfo,
    PathResult,
    PlanningResultInfo,
)
from raystar_interfaces.srv import GetRaystarPaths


def _field_names_and_types(message_type):
    return message_type.get_fields_and_field_types()


def test_action_schema_contract():
    expected_goal = {
        "map_id": "unique_identifier_msgs/UUID",
        "start": "geometry_msgs/PoseStamped",
        "goal": "geometry_msgs/PoseStamped",
        "search_mode": "uint8",
        "k": "int32",
        "max_path_length": "double",
        "allow_self_crossing": "boolean",
        "allow_unknown": "boolean",
        "include_debug": "boolean",
    }
    expected_result = {
        "success": "boolean",
        "result_info": "raystar_interfaces/PlanningResultInfo",
        "message": "string",
        "path_results": "sequence<raystar_interfaces/PathResult>",
        "debug_nodes": "sequence<raystar_interfaces/DebugNode>",
    }
    assert _field_names_and_types(PlanRaystarPaths.Goal) == expected_goal
    assert _field_names_and_types(PlanRaystarPaths.Result) == expected_result
    assert _field_names_and_types(PlanRaystarPaths.Feedback) == {}


def test_service_schema_contract():
    expected_request = {
        "map": "nav_msgs/OccupancyGrid",
        "start": "geometry_msgs/PoseStamped",
        "goal": "geometry_msgs/PoseStamped",
        "search_mode": "uint8",
        "k": "int32",
        "max_path_length": "double",
        "allow_self_crossing": "boolean",
        "allow_unknown": "boolean",
        "include_debug": "boolean",
    }
    expected_response = {
        "success": "boolean",
        "result_info": "raystar_interfaces/PlanningResultInfo",
        "message": "string",
        "path_results": "sequence<raystar_interfaces/PathResult>",
        "debug_nodes": "sequence<raystar_interfaces/DebugNode>",
    }
    assert _field_names_and_types(GetRaystarPaths.Request) == expected_request
    assert _field_names_and_types(GetRaystarPaths.Response) == expected_response


def test_goal_set_action_schema_contract():
    assert _field_names_and_types(PlanRaystarGoalSet.Goal) == {
        "map_id": "unique_identifier_msgs/UUID",
        "start": "geometry_msgs/PoseStamped",
        "goals": "sequence<geometry_msgs/PoseStamped>",
        "max_path_lengths": "sequence<double>",
        "allow_self_crossing": "boolean",
        "allow_unknown": "boolean",
        "include_debug": "boolean",
    }
    assert _field_names_and_types(PlanRaystarGoalSet.Result) == {
        "success": "boolean",
        "requested_start": "geometry_msgs/PoseStamped",
        "requested_allow_self_crossing": "boolean",
        "requested_allow_unknown": "boolean",
        "result_info": "raystar_interfaces/MultiGoalPlanningResultInfo",
        "message": "string",
        "goal_results": "sequence<raystar_interfaces/GoalPathResult>",
        "debug_nodes": "sequence<raystar_interfaces/DebugNode>",
    }
    assert _field_names_and_types(PlanRaystarGoalSet.Feedback) == {}


def test_transition_graph_action_schema_contract():
    assert _field_names_and_types(BuildRaystarTransitionGraph.Goal) == {
        "map_id": "unique_identifier_msgs/UUID",
        "expected_environment_id": "unique_identifier_msgs/UUID",
        "tether_configurations": "sequence<nav_msgs/Path>",
        "transition_pairs": "sequence<raystar_interfaces/ConfigurationTransitionPair>",
        "allow_unknown": "boolean",
    }
    assert _field_names_and_types(BuildRaystarTransitionGraph.Result) == {
        "success": "boolean",
        "map_id": "unique_identifier_msgs/UUID",
        "environment_id": "unique_identifier_msgs/UUID",
        "status": "uint8",
        "requested_transition_count": "uint32",
        "completed_transition_count": "uint32",
        "message": "string",
        "transitions": "sequence<raystar_interfaces/HomotopyTransitionResult>",
    }
    assert _field_names_and_types(BuildRaystarTransitionGraph.Feedback) == {
        "requested_transition_count": "uint32",
        "completed_transition_count": "uint32",
        "stage": "string",
    }


def test_message_schema_contract():
    assert _field_names_and_types(ConfigurationTransitionPair) == {
        "from_configuration": "uint32",
        "to_configuration": "uint32",
    }
    assert _field_names_and_types(HomotopyTransitionResult) == {
        "pair": "raystar_interfaces/ConfigurationTransitionPair",
        "status": "uint8",
        "path": "nav_msgs/Path",
        "path_length": "double",
        "collision_free": "boolean",
        "homotopy_preserved": "boolean",
        "locally_shortest": "boolean",
        "triangle_occurrences": "sequence<uint32>",
        "message": "string",
    }
    assert _field_names_and_types(DebugNode) == {
        "index": "int32",
        "g_cost": "double",
        "f_cost": "double",
    }
    assert _field_names_and_types(PathResult) == {
        "path": "nav_msgs/Path",
        "topology_path": "nav_msgs/Path",
        "cost": "double",
    }
    assert _field_names_and_types(MapStatus) == {
        "header": "std_msgs/Header",
        "map_id": "unique_identifier_msgs/UUID",
        "generation": "uint64",
        "width": "uint32",
        "height": "uint32",
        "resolution": "float",
        "occupied_threshold": "uint32",
        "environment_identity_version": "uint32",
        "occupancy_semantics_version": "uint32",
        "geometry_semantics_version": "uint32",
        "topology_semantics_version": "uint32",
    }
    assert _field_names_and_types(PlanningResultInfo) == {
        "map_id": "unique_identifier_msgs/UUID",
        "environment_id": "unique_identifier_msgs/UUID",
        "status": "uint8",
        "limits_reached": "uint16",
        "request_satisfied": "boolean",
        "search_complete": "boolean",
        "output_complete": "boolean",
        "debug_requested": "boolean",
        "debug_output_complete": "boolean",
        "requested_path_count": "uint32",
        "found_path_count": "uint32",
        "returned_path_count": "uint32",
        "expanded_nodes": "uint64",
        "map_time_ms": "double",
        "plan_time_ms": "double",
        "search_mode": "uint8",
        "requested_max_path_length": "double",
        "cost_bound_exhausted": "boolean",
    }
    assert _field_names_and_types(GoalPathResult) == {
        "goal_index": "uint32",
        "goal": "geometry_msgs/PoseStamped",
        "requested_max_path_length": "double",
        "success": "boolean",
        "result_info": "raystar_interfaces/PlanningResultInfo",
        "message": "string",
        "path_results": "sequence<raystar_interfaces/PathResult>",
    }
    assert _field_names_and_types(MultiGoalPlanningResultInfo) == {
        "map_id": "unique_identifier_msgs/UUID",
        "environment_id": "unique_identifier_msgs/UUID",
        "status": "uint8",
        "limits_reached": "uint16",
        "request_satisfied": "boolean",
        "search_complete": "boolean",
        "output_complete": "boolean",
        "debug_requested": "boolean",
        "debug_output_complete": "boolean",
        "requested_goal_count": "uint32",
        "returned_goal_count": "uint32",
        "completed_goal_count": "uint32",
        "goals_with_paths": "uint32",
        "found_path_count": "uint32",
        "returned_path_count": "uint32",
        "expanded_nodes": "uint64",
        "map_time_ms": "double",
        "plan_time_ms": "double",
    }


def test_status_constants_contract():
    assert BuildRaystarTransitionGraph.Result.STATUS_UNKNOWN == 0
    assert BuildRaystarTransitionGraph.Result.STATUS_COMPLETE == 1
    assert BuildRaystarTransitionGraph.Result.STATUS_INVALID_REQUEST == 2
    assert BuildRaystarTransitionGraph.Result.STATUS_CANCELLED == 3
    assert BuildRaystarTransitionGraph.Result.STATUS_FAILED == 4
    assert BuildRaystarTransitionGraph.Result.STATUS_TIMEOUT == 5
    assert HomotopyTransitionResult.STATUS_SUCCESS == 0
    assert HomotopyTransitionResult.STATUS_INVALID_REFERENCE == 1
    assert HomotopyTransitionResult.STATUS_NO_CORRIDOR == 2
    assert HomotopyTransitionResult.STATUS_STOPPED == 3
    assert HomotopyTransitionResult.STATUS_FAILURE == 4
    assert PlanningResultInfo.STATUS_UNKNOWN == 0
    assert PlanningResultInfo.STATUS_COMPLETE == 1
    assert PlanningResultInfo.STATUS_FEWER_PATHS == 2
    assert PlanningResultInfo.STATUS_NO_PATH == 3
    assert PlanningResultInfo.STATUS_PARTIAL_SEARCH == 4
    assert PlanningResultInfo.STATUS_PARTIAL_OUTPUT == 5
    assert PlanningResultInfo.STATUS_INVALID_REQUEST == 6
    assert PlanningResultInfo.STATUS_INVALID_CONFIGURATION == 7
    assert PlanningResultInfo.STATUS_BUSY == 8
    assert PlanningResultInfo.STATUS_CANCELLED == 9
    assert PlanningResultInfo.STATUS_FAILED == 10
    assert PlanningResultInfo.LIMIT_NONE == 0
    assert PlanningResultInfo.LIMIT_TIMEOUT == 1
    assert PlanningResultInfo.LIMIT_MAX_NODES == 2
    assert PlanningResultInfo.LIMIT_MAX_PATH_POINTS == 4
    assert PlanningResultInfo.LIMIT_MAX_RESPONSE_BYTES == 8
    assert PlanningResultInfo.LIMIT_CANCELLED == 16
    assert PlanningResultInfo.LIMIT_MAX_PATHS == 32
    assert PlanRaystarPaths.Goal.SEARCH_MODE_TOP_K == 0
    assert PlanRaystarPaths.Goal.SEARCH_MODE_ALL_WITHIN_LENGTH == 1
    assert GetRaystarPaths.Request.SEARCH_MODE_TOP_K == 0
    assert GetRaystarPaths.Request.SEARCH_MODE_ALL_WITHIN_LENGTH == 1
