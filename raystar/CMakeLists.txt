cmake_minimum_required(VERSION 3.8)
project(raystar)

# C++17
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# Find dependencies
find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(visualization_msgs REQUIRED)
find_package(raystar_interfaces REQUIRED)
find_package(CGAL REQUIRED)

# Core algorithm library (no ROS dependency in the algorithm itself)
add_library(raystar_core SHARED
  src/polymap.cpp
  src/raystar_core.cpp
)
target_include_directories(raystar_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_compile_features(raystar_core PUBLIC cxx_std_17)
target_link_libraries(raystar_core CGAL::CGAL)
# Avoid CGAL linkage issues
if(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE)
  set_property(TARGET raystar_core PROPERTY CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE ON)
endif()

# ROS2 node executable
add_executable(raystar_node
  src/raystar_node_main.cpp
  src/raystar_node.cpp
)
target_include_directories(raystar_node PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_link_libraries(raystar_node raystar_core)
ament_target_dependencies(raystar_node
  rclcpp
  nav_msgs
  geometry_msgs
  visualization_msgs
  raystar_interfaces
)

# Install
install(DIRECTORY include/
  DESTINATION include
)

install(TARGETS raystar_core
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

install(TARGETS raystar_node
  DESTINATION lib/${PROJECT_NAME}
)

# Tests
if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)

  ament_add_gtest(test_coordinate_utils
    test/test_coordinate_utils.cpp
  )
  target_include_directories(test_coordinate_utils PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )
  target_link_libraries(test_coordinate_utils raystar_core)

  ament_add_gtest(test_polymap
    test/test_polymap.cpp
  )
  target_include_directories(test_polymap PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )
  target_link_libraries(test_polymap raystar_core)

  ament_add_gtest(test_raystar_core
    test/test_raystar_core.cpp
  )
  target_include_directories(test_raystar_core PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )
  target_link_libraries(test_raystar_core raystar_core)

  ament_add_gtest(test_integration
    test/test_integration.cpp
  )
  target_include_directories(test_integration PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )
  target_compile_definitions(test_integration PRIVATE
    RAYSTAR_NODE_PATH="$<TARGET_FILE:raystar_node>"
  )
  target_link_libraries(test_integration raystar_core)
  ament_target_dependencies(test_integration
    rclcpp
    nav_msgs
    geometry_msgs
    visualization_msgs
    raystar_interfaces
  )
endif()

ament_package()
