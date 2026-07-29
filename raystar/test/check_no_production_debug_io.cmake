if(NOT DEFINED RAYSTAR_SOURCE_DIR)
  message(FATAL_ERROR "RAYSTAR_SOURCE_DIR is required")
endif()

file(GLOB_RECURSE production_files LIST_DIRECTORIES false
  "${RAYSTAR_SOURCE_DIR}/src/*.cpp"
  "${RAYSTAR_SOURCE_DIR}/src/*.h"
  "${RAYSTAR_SOURCE_DIR}/include/*.h"
  "${RAYSTAR_SOURCE_DIR}/include/*.hpp"
)

set(forbidden_debug_tokens
  "/tmp/"
  ".svg"
  "#include <fstream>"
  "std::fstream"
  "std::ofstream"
  "std::clog"
  "std::cout"
  "std::cerr"
  "fopen("
  "vis_call_count"
  "scope_call_count"
  "dumpSVG"
)

foreach(source_file IN LISTS production_files)
  file(READ "${source_file}" source_contents)
  foreach(forbidden_token IN LISTS forbidden_debug_tokens)
    string(FIND "${source_contents}" "${forbidden_token}" token_position)
    if(NOT token_position EQUAL -1)
      message(FATAL_ERROR
        "Production source file ${source_file} contains forbidden default debug I/O token: ${forbidden_token}")
    endif()
  endforeach()
endforeach()
