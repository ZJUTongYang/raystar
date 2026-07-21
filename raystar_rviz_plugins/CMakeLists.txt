cmake_minimum_required(VERSION 3.8)
project(raystar_rviz_plugins)

set(CMAKE_AUTOMOC ON)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(ament_cmake REQUIRED)
find_package(rviz_common REQUIRED)
find_package(rviz_rendering REQUIRED)
find_package(pluginlib REQUIRED)
find_package(rclcpp REQUIRED)
find_package(raystar_interfaces REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(Qt5 REQUIRED COMPONENTS Widgets)

set(rviz_plugins_deps
  rviz_common
  rviz_rendering
  pluginlib
  rclcpp
  raystar_interfaces
  nav_msgs
  geometry_msgs
)

add_library(raystar_rviz_plugins SHARED
  src/raystar_panel.cpp
  include/raystar_rviz_plugins/raystar_panel.h
)
target_include_directories(raystar_rviz_plugins PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
target_compile_features(raystar_rviz_plugins PUBLIC cxx_std_17)
ament_target_dependencies(raystar_rviz_plugins ${rviz_plugins_deps})
target_link_libraries(raystar_rviz_plugins Qt5::Widgets)

# Plugin description
pluginlib_export_plugin_description_file(rviz_common raystar_rviz_plugins.xml)

install(DIRECTORY include/
  DESTINATION include
)

install(TARGETS raystar_rviz_plugins
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

ament_package()
