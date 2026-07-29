cmake_minimum_required(VERSION 3.8)
project(raystar_relocated_install_smoke LANGUAGES CXX)

find_package(ament_cmake REQUIRED)

# This one find_package call must provide the complete public dependency set,
# including Qt5::Widgets, to an out-of-tree consumer.
find_package(raystar_rviz_plugins REQUIRED)
if(NOT TARGET raystar_rviz_plugins::raystar_rviz_plugins)
  message(FATAL_ERROR "raystar_rviz_plugins did not export its library target")
endif()
if(NOT TARGET Qt5::Widgets)
  message(FATAL_ERROR "raystar_rviz_plugins did not export its Qt Widgets dependency")
endif()

add_executable(export_target_smoke export_target_smoke.cpp)
target_compile_features(export_target_smoke PRIVATE cxx_std_17)
target_link_libraries(export_target_smoke PRIVATE
  raystar_rviz_plugins::raystar_rviz_plugins
)

# Deliberately do not link the Raystar plugin into this executable.  It must be
# discovered through the moved ament index and loaded from the plugin XML.
add_executable(plugin_loader_smoke plugin_loader_smoke.cpp)
target_compile_features(plugin_loader_smoke PRIVATE cxx_std_17)
ament_target_dependencies(plugin_loader_smoke
  pluginlib
  rviz_common
)
# ament_target_dependencies() uses CMake's plain target_link_libraries()
# signature on Humble, so this call must use the same signature.
target_link_libraries(plugin_loader_smoke Qt5::Widgets)
