foreach(_required
    RAYSTAR_PROFILE_BINARY
    RAYSTAR_PROFILE_SOURCE
    RAYSTAR_CMAKE_FILE
    RAYSTAR_SOURCE_TESTMAP)
  if(NOT DEFINED ${_required} OR "${${_required}}" STREQUAL "")
    message(FATAL_ERROR "${_required} must identify a profiling contract input")
  endif()
endforeach()

foreach(_required_file
    "${RAYSTAR_PROFILE_BINARY}"
    "${RAYSTAR_PROFILE_SOURCE}"
    "${RAYSTAR_CMAKE_FILE}"
    "${RAYSTAR_SOURCE_TESTMAP}")
  if(NOT EXISTS "${_required_file}")
    message(FATAL_ERROR "Profiling contract input does not exist: ${_required_file}")
  endif()
endforeach()

file(READ "${RAYSTAR_CMAKE_FILE}" _cmake_text)
file(READ "${RAYSTAR_PROFILE_SOURCE}" _source_text)

foreach(_forbidden RAYSTAR_PROFILE_TESTMAP_PATH "/proc/self/exe")
  string(FIND "${_cmake_text}" "${_forbidden}" _cmake_match)
  string(FIND "${_source_text}" "${_forbidden}" _source_match)
  if(NOT _cmake_match EQUAL -1 OR NOT _source_match EQUAL -1)
    message(FATAL_ERROR
      "Profiling source/build configuration retains forbidden fallback '${_forbidden}'")
  endif()
endforeach()

foreach(_required_source_contract
    "ament_index_cpp::get_package_share_directory(\"raystar\")"
    "argument == \"--testmap\""
    "options.testmap_path = require_value(argument)")
  string(FIND "${_source_text}" "${_required_source_contract}" _source_match)
  if(_source_match EQUAL -1)
    message(FATAL_ERROR
      "Profiling source is missing contract '${_required_source_contract}'")
  endif()
endforeach()

# file(STRINGS) extracts printable binary payload without requiring a
# platform-specific `strings` executable.  The installed-relative suffix is
# intentionally present; only the configure-time source-tree path is banned.
file(STRINGS "${RAYSTAR_PROFILE_BINARY}" _binary_strings)
foreach(_binary_string IN LISTS _binary_strings)
  string(FIND "${_binary_string}" "${RAYSTAR_SOURCE_TESTMAP}" _source_map_match)
  if(NOT _source_map_match EQUAL -1)
    message(FATAL_ERROR
      "raystar_profile embeds source-tree map path '${RAYSTAR_SOURCE_TESTMAP}'")
  endif()
  foreach(_forbidden RAYSTAR_PROFILE_TESTMAP_PATH "/proc/self/exe")
    string(FIND "${_binary_string}" "${_forbidden}" _binary_match)
    if(NOT _binary_match EQUAL -1)
      message(FATAL_ERROR
        "raystar_profile embeds forbidden fallback '${_forbidden}'")
    endif()
  endforeach()
endforeach()
