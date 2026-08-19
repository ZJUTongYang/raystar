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
find_package(ament_index_cpp REQUIRED)
find_package(rclcpp REQUIRED)
find_package(rclcpp_action REQUIRED)
find_package(rcl_interfaces REQUIRED)
find_package(nav_msgs REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(std_msgs REQUIRED)
find_package(visualization_msgs REQUIRED)
find_package(raystar_interfaces REQUIRED)
find_package(unique_identifier_msgs REQUIRED)
find_package(Boost REQUIRED)
find_package(CGAL REQUIRED)

# CGAL 5.x's ConfigVersion file rejects a requested 5.x version when CGAL 6.x
# is installed.  Compare the discovered version ourselves so the lower bound
# remains valid across CGAL major releases.
if(NOT DEFINED CGAL_VERSION OR CGAL_VERSION VERSION_LESS "5.4")
  message(FATAL_ERROR
    "Raystar requires CGAL >= 5.4; found '${CGAL_VERSION}'")
endif()
# The exact-constructions kernel used by Raystar requires CGAL's GMP/MPFR
# backend.  Keep this explicit even though CGAL normally discovers both while
# defining CGAL::CGAL.
if(CGAL_DISABLE_GMP OR NOT GMP_FOUND OR NOT MPFR_FOUND)
  message(FATAL_ERROR
    "Raystar requires CGAL with both GMP and MPFR support enabled")
endif()

# Package-private core algorithm library (no ROS dependency in the algorithm itself).
# It is linked statically into the node/tests and is intentionally not installed
# or exported as a downstream CMake API.
add_library(raystar_core STATIC
  src/polymap.cpp
  src/polymap_visibility.cpp
  src/polymap_shortening.cpp
  src/raystar_core.cpp
  src/raystar_core_planning.cpp
  src/visibility_validation.cpp
)
target_include_directories(raystar_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_compile_features(raystar_core PUBLIC cxx_std_17)
# PUBLIC only propagates build-tree requirements to this package's node/tests;
# raystar_core itself is not part of the installed package interface.
target_link_libraries(raystar_core PUBLIC CGAL::CGAL Boost::boost)
# Avoid CGAL linkage issues
if(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE)
  set_property(TARGET raystar_core PROPERTY CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE ON)
endif()

option(RAYSTAR_BUILD_PROFILING
  "Build the deterministic Raystar profiling and acceptance runner" OFF)
if(RAYSTAR_BUILD_PROFILING)
  add_executable(raystar_profile
    tools/raystar_profile.cpp
  )
  target_include_directories(raystar_profile PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/test
  )
  target_compile_features(raystar_profile PRIVATE cxx_std_17)
  target_link_libraries(raystar_profile raystar_core)
  ament_target_dependencies(raystar_profile
    ament_index_cpp
    geometry_msgs
    nav_msgs
  )
  install(TARGETS raystar_profile
    DESTINATION lib/${PROJECT_NAME}
  )
  install(PROGRAMS tools/summarize_profile.py
    DESTINATION lib/${PROJECT_NAME}
    RENAME raystar_profile_summary
  )
endif()

# ROS2 node executable
add_executable(raystar_node
  src/raystar_node_main.cpp
  src/raystar_node.cpp
  src/raystar_node_actions.cpp
  src/raystar_node_planning.cpp
  src/raystar_node_visualization.cpp
)
target_include_directories(raystar_node PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)
target_link_libraries(raystar_node raystar_core)
ament_target_dependencies(raystar_node
  rclcpp
  rclcpp_action
  rcl_interfaces
  nav_msgs
  geometry_msgs
  std_msgs
  visualization_msgs
  raystar_interfaces
  unique_identifier_msgs
)

# Install only the ROS node and runtime resources. The algorithm archive and
# internal headers remain in the build/source tree.
install(TARGETS raystar_node
  DESTINATION lib/${PROJECT_NAME}
)

install(DIRECTORY launch rviz
  DESTINATION share/${PROJECT_NAME}
  PATTERN "__pycache__" EXCLUDE
  PATTERN "*.pyc" EXCLUDE
)

install(DIRECTORY test/
  DESTINATION share/${PROJECT_NAME}/test
  FILES_MATCHING
    PATTERN "*.yaml"
    PATTERN "*.pgm"
    PATTERN "__pycache__" EXCLUDE
    PATTERN "*.pyc" EXCLUDE
)

# Tests
if(BUILD_TESTING)
  find_package(ament_cmake_gtest REQUIRED)
  find_package(ament_lint_auto REQUIRED)
  # Humble and Jazzy install different optional flake8 plugins.  Pin the
  # package policy so identical Python sources receive identical checks in
  # every supported ROS distribution.
  set(ament_cmake_flake8_CONFIG_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}/test/ament_flake8.ini")
  # Let the package's declared lint test dependencies register their checks
  # with ament_package().  Keeping this in BUILD_TESTING avoids pulling lint
  # tools into a production-only configure.
  ament_lint_auto_find_test_dependencies()

  ament_add_gtest(test_coordinate_utils
    test/test_coordinate_utils.cpp
  )
  target_include_directories(test_coordinate_utils PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )
  target_link_libraries(test_coordinate_utils raystar_core)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(test_coordinate_utils PRIVATE -frounding-math)
  endif()

  ament_add_gtest(test_conservative_path_length
    test/test_conservative_path_length.cpp
  )
  target_include_directories(test_conservative_path_length PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )
  target_link_libraries(test_conservative_path_length Boost::boost)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(test_conservative_path_length PRIVATE -frounding-math)
  endif()

  ament_add_gtest(test_published_path_order
    test/test_published_path_order.cpp
  )
  target_include_directories(test_published_path_order PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )
  ament_target_dependencies(test_published_path_order
    raystar_interfaces
  )

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
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )
  target_link_libraries(test_raystar_core raystar_core)

  ament_add_gtest(test_visibility_validation
    test/test_visibility_validation.cpp
  )
  target_include_directories(test_visibility_validation PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )
  target_link_libraries(test_visibility_validation raystar_core)

  ament_add_gtest(test_visibility_reference
    test/test_visibility_reference.cpp
  )
  target_include_directories(test_visibility_reference PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )
  target_link_libraries(test_visibility_reference raystar_core)

  ament_add_gtest(test_topov_regression
    test/test_topov_regression.cpp
  )
  target_include_directories(test_topov_regression PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )
  target_compile_definitions(test_topov_regression PRIVATE
    RAYSTAR_TESTMAP_PATH="${CMAKE_CURRENT_SOURCE_DIR}/test/testmap.pgm"
  )
  target_link_libraries(test_topov_regression raystar_core)

  # The public production headers keep planner state encapsulated.  A small
  # set of existing white-box regression tests still probes invariants and
  # failure injection points directly; expose those members only in these
  # test translation units.
  foreach(_raystar_white_box_test
      test_polymap
      test_raystar_core
      test_visibility_validation
      test_visibility_reference
      test_topov_regression)
    target_compile_definitions(${_raystar_white_box_test} PRIVATE
      RAYSTAR_TESTING=1
    )
  endforeach()

  # Independent final-path obstacle-collision and self-intersection oracle.
  # This target intentionally does not link raystar_core so its geometry cannot
  # accidentally reuse production coordinate, visibility, or pruning helpers.
  ament_add_gtest(test_path_collision_oracle
    test/test_path_collision_oracle.cpp
  )
  ament_target_dependencies(test_path_collision_oracle
    geometry_msgs
    nav_msgs
  )

  # Fixed-seed small-grid generator and independent four-neighbor BFS oracle.
  # Like the path oracle, this target deliberately avoids raystar_core so a
  # production reachability regression cannot be reproduced in the reference.
  ament_add_gtest(test_small_grid_property
    test/test_small_grid_property.cpp
  )
  ament_target_dependencies(test_small_grid_property
    nav_msgs
  )

  ament_add_gtest(test_map_identity
    test/test_map_identity.cpp
  )
  ament_target_dependencies(test_map_identity
    nav_msgs
    raystar_interfaces
    unique_identifier_msgs
  )

  ament_add_gtest(test_transition_environment_cache
    test/test_transition_environment_cache.cpp
  )
  target_include_directories(test_transition_environment_cache PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  )
  target_link_libraries(test_transition_environment_cache raystar_core)

  ament_add_gtest(test_transition_progress_reporter
    test/test_transition_progress_reporter.cpp
  )
  target_include_directories(test_transition_progress_reporter PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )

  add_test(NAME test_no_production_debug_io
    COMMAND ${CMAKE_COMMAND}
      -DRAYSTAR_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_no_production_debug_io.cmake
  )
  set_tests_properties(test_no_production_debug_io PROPERTIES LABELS "lint")

  add_test(NAME test_rviz_config_contract
    COMMAND ${CMAKE_COMMAND}
      -DRAYSTAR_RVIZ_CONFIG=${CMAKE_CURRENT_SOURCE_DIR}/rviz/raystar_test.rviz
      -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_rviz_config.cmake
  )
  set_tests_properties(test_rviz_config_contract PROPERTIES LABELS "lint;rviz")

  if(RAYSTAR_BUILD_PROFILING)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    add_test(NAME test_profile_relocation_contract
      COMMAND ${CMAKE_COMMAND}
        -DRAYSTAR_PROFILE_BINARY=$<TARGET_FILE:raystar_profile>
        -DRAYSTAR_PROFILE_SOURCE=${CMAKE_CURRENT_SOURCE_DIR}/tools/raystar_profile.cpp
        -DRAYSTAR_CMAKE_FILE=${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
        -DRAYSTAR_SOURCE_TESTMAP=${CMAKE_CURRENT_SOURCE_DIR}/test/testmap.pgm
        -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_profile_relocation.cmake
    )
    set_tests_properties(test_profile_relocation_contract PROPERTIES
      LABELS "profiling;packaging"
      TIMEOUT 30
    )

    add_test(NAME test_profile_smoke
      COMMAND $<TARGET_FILE:raystar_profile>
        --scenario open_256
        --k-values 1
        --warmups 0
        --iterations 1
    )
    set_tests_properties(test_profile_smoke PROPERTIES
      LABELS "profiling"
      TIMEOUT 30
    )

    add_test(NAME test_profile_all_within_length_smoke
      COMMAND $<TARGET_FILE:raystar_profile>
        --scenario open_256
        --mode all-within-length
        --max-path-length 10
        --warmups 0
        --iterations 1
    )
    set_tests_properties(test_profile_all_within_length_smoke PROPERTIES
      LABELS "profiling"
      TIMEOUT 30
    )

    add_test(NAME test_profile_multi_goal_smoke
      COMMAND $<TARGET_FILE:raystar_profile>
        --scenario open_256
        --mode multi-goal
        --max-path-length 10
        --goal-count 3
        --warmups 0
        --iterations 1
    )
    set_tests_properties(test_profile_multi_goal_smoke PROPERTIES
      LABELS "profiling"
      TIMEOUT 30
    )

    add_test(NAME test_profile_bundled_map_lookup
      COMMAND $<TARGET_FILE:raystar_profile>
        --scenario bundled_testmap_50
        --k-values 1
        --warmups 0
        --iterations 1
    )
    set_tests_properties(test_profile_bundled_map_lookup PROPERTIES
      LABELS "profiling;packaging"
      TIMEOUT 30
    )

    add_test(NAME test_profile_explicit_map_override
      COMMAND ${CMAKE_COMMAND} -E env AMENT_PREFIX_PATH=
        $<TARGET_FILE:raystar_profile>
        --scenario bundled_testmap_50
        --testmap ${CMAKE_CURRENT_SOURCE_DIR}/test/testmap.pgm
        --k-values 1
        --warmups 0
        --iterations 1
    )
    set_tests_properties(test_profile_explicit_map_override PROPERTIES
      LABELS "profiling;packaging"
      TIMEOUT 30
    )

    function(add_profile_failure_test test_name)
      add_test(NAME ${test_name}
        COMMAND $<TARGET_FILE:raystar_profile> ${ARGN}
      )
      set_tests_properties(${test_name} PROPERTIES
        LABELS "profiling"
        TIMEOUT 30
        WILL_FAIL TRUE
      )
    endfunction()

    add_profile_failure_test(test_profile_reject_missing_bound
      --scenario open_256
      --mode all-within-length
      --warmups 0
      --iterations 1
    )
    add_profile_failure_test(test_profile_reject_nan_bound
      --scenario open_256
      --mode all-within-length
      --max-path-length nan
      --warmups 0
      --iterations 1
    )
    add_profile_failure_test(test_profile_reject_bounded_k_values
      --scenario open_256
      --mode all-within-length
      --max-path-length 10
      --k-values 1
      --warmups 0
      --iterations 1
    )
    add_profile_failure_test(test_profile_reject_invalid_goal_count
      --scenario open_256
      --mode multi-goal
      --max-path-length 10
      --goal-count 1
      --warmups 0
      --iterations 1
    )
    add_profile_failure_test(test_profile_reject_path_point_limit
      --scenario open_256
      --max-path-points 1
      --warmups 0
      --iterations 1
    )
    add_profile_failure_test(test_profile_reject_multi_goal_admission
      --scenario open_256
      --mode multi-goal
      --max-path-length 10
      --goal-count 3
      --max-multi-goal-count 2
      --warmups 0
      --iterations 1
    )
    add_profile_failure_test(test_profile_immediate_timeout
      --scenario open_256
      --mode all-within-length
      --max-path-length 400
      --timeout-ms 0
      --warmups 0
      --iterations 1
    )

    add_test(NAME test_profile_summary
      COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test/test_profile_summary.py
        ${CMAKE_CURRENT_SOURCE_DIR}/tools/summarize_profile.py
    )
    set_tests_properties(test_profile_summary PROPERTIES
      LABELS "profiling"
      TIMEOUT 30
    )
  endif()

  ament_add_gtest(test_integration
    test/test_integration.cpp
  )
  target_include_directories(test_integration PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  )
  target_link_libraries(test_integration Boost::boost)
  target_compile_definitions(test_integration PRIVATE
    RAYSTAR_NODE_PATH="$<TARGET_FILE:raystar_node>"
  )
  ament_target_dependencies(test_integration
    rclcpp
    rclcpp_action
    rcl_interfaces
    nav_msgs
    geometry_msgs
    visualization_msgs
    raystar_interfaces
  )
endif()

ament_package()
