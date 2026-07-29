"""Contract tests for the wire-facing Raystar ROS interfaces.

The node and RViz plugin are intentionally built against these generated
types.  This test keeps an accidental field rename, reorder, or type change
from silently becoming a wire-compatibility break.  The checks use the
generated Python metadata rather than duplicating the IDL parser, so they also
verify that rosidl generated the expected public schema.
"""

from raystar_interfaces.action import PlanRaystarPaths
from raystar_interfaces.msg import DebugNode, MapStatus, PathResult, PlanningResultInfo
from raystar_interfaces.srv import GetRaystarPaths


def _field_names_and_types(message_type):
    return message_type.get_fields_and_field_types()


def test_action_schema_contract():
    expected_goal = {
        "map_id": "unique_identifier_msgs/UUID",
        "start": "geometry_msgs/PoseStamped",
        "goal": "geometry_msgs/PoseStamped",
        "k": "int32",
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
        "k": "int32",
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


def test_message_schema_contract():
    assert _field_names_and_types(DebugNode) == {
        "index": "int32",
        "g_cost": "double",
        "f_cost": "double",
    }
    assert _field_names_and_types(PathResult) == {
        "path": "nav_msgs/Path",
        "cost": "double",
    }
    assert _field_names_and_types(MapStatus) == {
        "header": "std_msgs/Header",
        "map_id": "unique_identifier_msgs/UUID",
        "generation": "uint64",
        "width": "uint32",
        "height": "uint32",
        "resolution": "float",
    }
    assert _field_names_and_types(PlanningResultInfo) == {
        "map_id": "unique_identifier_msgs/UUID",
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
    }


def test_status_constants_contract():
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
