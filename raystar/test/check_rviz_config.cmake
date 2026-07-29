if(NOT DEFINED RAYSTAR_RVIZ_CONFIG)
  message(FATAL_ERROR "RAYSTAR_RVIZ_CONFIG is required")
endif()

file(READ "${RAYSTAR_RVIZ_CONFIG}" content)
string(
  REGEX MATCHALL
  "Class:[ \t]+raystar_rviz_plugins[^ \t\r\n]*"
  panel_entries
  "${content}"
)
list(LENGTH panel_entries panel_count)

if(NOT panel_count EQUAL 1)
  message(FATAL_ERROR
    "Expected exactly one Raystar Panel in the bundled RViz config, got "
    "${panel_count}: ${panel_entries}")
endif()

if(NOT content MATCHES
    "Class:[ \t]+raystar_rviz_plugins/RaystarPanel([ \t\r\n]|$)")
  message(FATAL_ERROR
    "Bundled RViz config must use raystar_rviz_plugins/RaystarPanel")
endif()

if(content MATCHES "Class:[ \t]+raystar_rviz_plugins::RaystarPanel")
  message(FATAL_ERROR
    "Bundled RViz config contains the obsolete C++ type lookup ID")
endif()

foreach(required_setting
    "action_name:[ \t]+/raystar/plan_paths"
    "allow_self_crossing:[ \t]+false"
    "allow_unknown:[ \t]+false")
  if(NOT content MATCHES "${required_setting}")
    message(FATAL_ERROR
      "Bundled RViz config is missing required setting: ${required_setting}")
  endif()
endforeach()
