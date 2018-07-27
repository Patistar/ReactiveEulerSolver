set(CMAKE_REQUIRED_FLAGS "-std=c++${FUB_CORE_CXX_STANDARD}")
check_cxx_source_compiles("
#include <variant>
int main() { std::variant<int> _; }
"
        FUB_WITH_STD_VARIANT)

if (NOT FUB_WITH_STD_VARIANT)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/variant/.git)
    execute_process(COMMAND git submodule update --init --remote --depth=1 -- third-party/variant
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  add_library(mpark_variant INTERFACE)
  target_include_directories(mpark_variant INTERFACE
          ${CMAKE_CURRENT_SOURCE_DIR}/third-party/variant/include)
endif()
