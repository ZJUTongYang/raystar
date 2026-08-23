cmake_minimum_required(VERSION 3.8)
project(raystar_interfaces)

find_package(ament_cmake REQUIRED)
find_package(rosidl_default_generators REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(std_msgs REQUIRED)
find_package(unique_identifier_msgs REQUIRED)
find_package(action_msgs REQUIRED)

rosidl_generate_interfaces(${PROJECT_NAME}
  "action/BuildRaystarTransitionGraph.action"
  "action/PlanRaystarGoalSet.action"
  "action/PlanRaystarPaths.action"
  "msg/ReferenceTransitionPair.msg"
  "msg/DebugNode.msg"
  "msg/GoalPathResult.msg"
  "msg/HomotopyTransitionResult.msg"
  "msg/MapStatus.msg"
  "msg/MultiGoalPlanningResultInfo.msg"
  "msg/PathResult.msg"
  "msg/PlanningResultInfo.msg"
  "srv/GetRaystarPaths.srv"
  DEPENDENCIES nav_msgs geometry_msgs std_msgs unique_identifier_msgs action_msgs
)

install(DIRECTORY include/
  DESTINATION include
)
# rosidl exports <prefix>/include/${PROJECT_NAME} as the modern target's
# include root.  Install the helper there as well so consumers can use the
# conventional <raystar_interfaces/map_identity.hpp> spelling.  The directory
# install above preserves the original <map_identity.hpp> path for source and
# binary compatibility with this development branch.
install(FILES
  include/raystar_interfaces/environment_identity.hpp
  include/raystar_interfaces/map_identity.hpp
  DESTINATION include/${PROJECT_NAME}/${PROJECT_NAME}
)
ament_export_include_directories(include)
ament_export_dependencies(
  nav_msgs
  geometry_msgs
  std_msgs
  unique_identifier_msgs
  action_msgs
)

if(BUILD_TESTING)
  find_package(ament_cmake_pytest REQUIRED)
  find_package(ament_lint_auto REQUIRED)
  ament_lint_auto_find_test_dependencies()
  ament_add_pytest_test(
    test_interface_schema
    test/test_interface_schema.py
  )
endif()

ament_package()
