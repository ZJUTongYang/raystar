"""Launch the bundled Ray* demo stack.

The launch file is deliberately useful both as a one-command demo and as a
small bringup template.  It starts a lifecycle-managed Nav2 map server, the
Raystar node, and (optionally) RViz.  All large planner/resource settings are
forwarded as node parameters, while ``map_topic`` and ``action_name`` are
normal ROS remappings so a namespaced deployment can choose its own endpoints.

For a headless server, use ``start_rviz:=false``.  For a deployment that
already owns a map server, use ``start_map_server:=false`` and point
``map_topic`` at that server's OccupancyGrid topic.  The default RViz config
uses ``/raystar/plan_paths``; when a namespace or custom action endpoint is
selected, update the panel's persisted Action field accordingly.
"""

from launch import LaunchDescription
from launch.actions import DeclareLaunchArgument
from launch.conditions import IfCondition
from launch.substitutions import LaunchConfiguration, PathJoinSubstitution
from launch_ros.actions import Node
from launch_ros.parameter_descriptions import ParameterValue
from launch_ros.substitutions import FindPackageShare


def _bool_parameter(name):
    """Return a typed launch substitution for a boolean ROS parameter."""
    return ParameterValue(LaunchConfiguration(name), value_type=bool)


def _int_parameter(name):
    """Return a typed launch substitution for an integer ROS parameter."""
    return ParameterValue(LaunchConfiguration(name), value_type=int)


def generate_launch_description():
    namespace = LaunchConfiguration("namespace")
    map_topic = LaunchConfiguration("map_topic")
    action_name = LaunchConfiguration("action_name")
    goal_set_action_name = LaunchConfiguration("goal_set_action_name")
    transition_action_name = LaunchConfiguration("transition_action_name")
    use_sim_time = _bool_parameter("use_sim_time")

    default_map = PathJoinSubstitution(
        [FindPackageShare("raystar"), "test", "testmap.yaml"]
    )
    default_rviz = PathJoinSubstitution(
        [FindPackageShare("raystar"), "rviz", "raystar_test.rviz"]
    )

    arguments = [
        DeclareLaunchArgument(
            "map_yaml",
            default_value=default_map,
            description="Map YAML consumed by nav2_map_server.",
        ),
        DeclareLaunchArgument(
            "namespace",
            default_value="",
            description="Namespace for map server, lifecycle manager, and Raystar.",
        ),
        DeclareLaunchArgument(
            "map_topic",
            default_value="/map",
            description="OccupancyGrid topic cached by Raystar and published by the map server.",
        ),
        DeclareLaunchArgument(
            "action_name",
            default_value="~/plan_paths",
            description=(
                "Action endpoint remapped from Raystar's private ~/plan_paths name. "
                "Use an absolute name for a shared endpoint."
            ),
        ),
        DeclareLaunchArgument(
            "goal_set_action_name",
            default_value="~/plan_goal_set",
            description=(
                "Action endpoint remapped from Raystar's private ~/plan_goal_set name."
            ),
        ),
        DeclareLaunchArgument(
            "transition_action_name",
            default_value="~/build_transition_graph",
            description=(
                "Action endpoint remapped from Raystar's private "
                "~/build_transition_graph name."
            ),
        ),
        DeclareLaunchArgument(
            "rviz_config",
            default_value=default_rviz,
            description="RViz2 configuration file.",
        ),
        DeclareLaunchArgument(
            "start_map_server",
            default_value="true",
            description="Start and activate the bundled Nav2 map server.",
        ),
        DeclareLaunchArgument(
            "start_rviz",
            default_value="true",
            description="Start RViz2 with rviz_config.",
        ),
        DeclareLaunchArgument(
            "use_sim_time",
            default_value="false",
            description="Use the /clock topic for all launched nodes.",
        ),
        DeclareLaunchArgument(
            "enable_legacy_map_service",
            default_value="false",
            description="Expose the compatibility full-map Service (larger DDS samples).",
        ),
        DeclareLaunchArgument("occupied_threshold", default_value="99"),
        DeclareLaunchArgument("max_k", default_value="100"),
        DeclareLaunchArgument(
            "max_cost_bounded_paths",
            default_value="1000",
            description="Per-goal path cap for exhaustive bounded enumeration.",
        ),
        DeclareLaunchArgument(
            "max_multi_goal_count",
            default_value="128",
            description="Goal-array admission limit for one shared-tree request.",
        ),
        DeclareLaunchArgument(
            "max_transition_configurations",
            default_value="4096",
            description="Tether-configuration admission limit for one UPS batch.",
        ),
        DeclareLaunchArgument(
            "max_transition_pairs",
            default_value="1000",
            description="Directed-pair admission limit for one UPS batch.",
        ),
        DeclareLaunchArgument("max_nodes", default_value="10000"),
        DeclareLaunchArgument(
            "planning_timeout_ms",
            default_value="5000",
            description="Cooperative deadline for each planning or UPS request.",
        ),
        DeclareLaunchArgument("max_map_cells", default_value="8388608"),
        DeclareLaunchArgument("max_map_bytes", default_value="536870912"),
        DeclareLaunchArgument("max_path_points", default_value="100000"),
        DeclareLaunchArgument("max_debug_nodes", default_value="0"),
        DeclareLaunchArgument("max_response_bytes", default_value="67108864"),
        DeclareLaunchArgument(
            "path_visualization_republish_period_ms",
            default_value="2000",
            description="Periodic cached path MarkerArray refresh; set 0 to disable.",
        ),
    ]

    map_server = Node(
        package="nav2_map_server",
        executable="map_server",
        name="map_server",
        namespace=namespace,
        output="screen",
        condition=IfCondition(LaunchConfiguration("start_map_server")),
        parameters=[
            {
                "yaml_filename": LaunchConfiguration("map_yaml"),
                "use_sim_time": use_sim_time,
            }
        ],
        remappings=[("map", map_topic)],
    )

    lifecycle_manager = Node(
        package="nav2_lifecycle_manager",
        executable="lifecycle_manager",
        name="lifecycle_manager_map",
        namespace=namespace,
        output="screen",
        condition=IfCondition(LaunchConfiguration("start_map_server")),
        parameters=[
            {
                "autostart": True,
                "use_sim_time": use_sim_time,
                "node_names": ["map_server"],
            }
        ],
    )

    raystar_node = Node(
        package="raystar",
        executable="raystar_node",
        name="raystar",
        namespace=namespace,
        output="screen",
        parameters=[
            {
                "map_topic": map_topic,
                "enable_legacy_map_service": _bool_parameter(
                    "enable_legacy_map_service"
                ),
                "occupied_threshold": _int_parameter("occupied_threshold"),
                "max_k": _int_parameter("max_k"),
                "max_cost_bounded_paths": _int_parameter(
                    "max_cost_bounded_paths"
                ),
                "max_multi_goal_count": _int_parameter("max_multi_goal_count"),
                "max_transition_configurations": _int_parameter(
                    "max_transition_configurations"
                ),
                "max_transition_pairs": _int_parameter("max_transition_pairs"),
                "max_nodes": _int_parameter("max_nodes"),
                "planning_timeout_ms": _int_parameter("planning_timeout_ms"),
                "max_map_cells": _int_parameter("max_map_cells"),
                "max_map_bytes": _int_parameter("max_map_bytes"),
                "max_path_points": _int_parameter("max_path_points"),
                "max_debug_nodes": _int_parameter("max_debug_nodes"),
                "max_response_bytes": _int_parameter("max_response_bytes"),
                "path_visualization_republish_period_ms": _int_parameter(
                    "path_visualization_republish_period_ms"
                ),
            }
        ],
        remappings=[
            ("~/plan_paths", action_name),
            ("~/plan_goal_set", goal_set_action_name),
            ("~/build_transition_graph", transition_action_name),
        ],
    )

    rviz = Node(
        package="rviz2",
        executable="rviz2",
        name="rviz2",
        output="screen",
        condition=IfCondition(LaunchConfiguration("start_rviz")),
        arguments=["-d", LaunchConfiguration("rviz_config")],
        parameters=[{"use_sim_time": use_sim_time}],
    )

    return LaunchDescription(arguments + [map_server, lifecycle_manager, raystar_node, rviz])
