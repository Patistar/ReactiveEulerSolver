# This script checks whether an implementation for std::optional is available.

include(CheckCXXSourceCompiles)

get_filename_component(_mod_dir ${CMAKE_CURRENT_LIST_FILE} PATH)

set(CMAKE_REQUIRED_FLAGS "-std=c++${FUB_CORE_CXX_STANDARD}")
check_cxx_source_compiles("
#include <optional>
int main() { std::optional<int> _; }
"
        FUB_CORE_USE_STD_OPTIONAL)

add_library(ReactiveEulerSolver_core_optional INTERFACE)
if (FUB_CORE_USE_STD_OPTIONAL)
  target_compile_features(ReactiveEulerSolver_core_optional INTERFACE cxx_std_17)
else()
  if (NOT EXISTS ${CMAKE_SOURCE_DIR}/core/third-party/optional-lite/.git)
    execute_process(COMMAND git submodule update --remote --init --depth=1 -- third-party/optional-lite
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/core)
  endif()
  add_subdirectory(${CMAKE_SOURCE_DIR}/core/third-party/optional-lite)
  target_link_libraries(ReactiveEulerSolver_core_optional INTERFACE optional-lite)
endif()

unset(_mod_dir)
