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
find_package(rviz_common 11.2 REQUIRED)
find_package(pluginlib 5.1 REQUIRED)
find_package(rclcpp REQUIRED)
find_package(rclcpp_action REQUIRED)
find_package(raystar_interfaces REQUIRED)
find_package(unique_identifier_msgs REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(Qt5 5.15 REQUIRED COMPONENTS Widgets)

set(rviz_plugins_deps
  rviz_common
  pluginlib
  rclcpp
  rclcpp_action
  raystar_interfaces
  unique_identifier_msgs
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
  EXPORT export_${PROJECT_NAME}
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)
  find_package(ament_lint_auto REQUIRED)
  ament_lint_auto_find_test_dependencies()

  # Some ROS 2 logging backends open a file as soon as a Qt/RViz test creates
  # its first node.  Keep the test hermetic when colcon is invoked without a
  # user-provided ROS_LOG_DIR (for example, on a read-only CI home).
  set(_raystar_plugin_test_log_dir
    "${CMAKE_CURRENT_BINARY_DIR}/test_ros_logs")
  file(MAKE_DIRECTORY "${_raystar_plugin_test_log_dir}")

  pluginlib_enable_plugin_testing(
    PACKAGE_NAME ${PROJECT_NAME}
    PLUGIN_CATEGORY rviz_common
    PLUGIN_DESCRIPTIONS raystar_rviz_plugins.xml
    PLUGIN_LIBRARIES raystar_rviz_plugins
    CMAKE_TARGET_VAR plugin_test_mock_target
    AMENT_PREFIX_PATH_VAR plugin_test_prefix
  )

  ament_add_gtest(test_raystar_panel_plugin
    test/test_raystar_panel_plugin.cpp
    SKIP_LINKING_MAIN_LIBRARIES
    ENV "QT_QPA_PLATFORM=offscreen"
        "ROS_LOG_DIR=${_raystar_plugin_test_log_dir}"
    APPEND_ENV "AMENT_PREFIX_PATH=${plugin_test_prefix}"
  )
  if(TARGET test_raystar_panel_plugin)
    ament_target_dependencies(test_raystar_panel_plugin
      pluginlib
      rviz_common
      rclcpp
      rclcpp_action
      raystar_interfaces
      unique_identifier_msgs
      nav_msgs
      geometry_msgs
    )
    target_include_directories(test_raystar_panel_plugin PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_link_libraries(test_raystar_panel_plugin
      raystar_rviz_plugins
      Qt5::Widgets
    )
    add_dependencies(test_raystar_panel_plugin
      raystar_rviz_plugins
      "${plugin_test_mock_target}"
    )
  endif()
endif()

# The panel header is installed as a public plugin API.  Export the target,
# include path, and all transitive ROS/Qt dependencies so a downstream CMake
# package can link to it without relying on this source tree's build paths.
ament_export_targets(export_${PROJECT_NAME} HAS_LIBRARY_TARGET)
ament_export_include_directories(include)
ament_export_libraries(raystar_rviz_plugins)
ament_export_dependencies(${rviz_plugins_deps})

ament_package(CONFIG_EXTRAS cmake/raystar_rviz_plugins-extras.cmake)
