# This script checks whether an implementation for std::variant is available.

get_filename_component(_mod_dir ${CMAKE_CURRENT_LIST_FILE} PATH)
message(STATUS "Checking for std::variant...")
try_compile(_has_std_variant
  ${CMAKE_BINARY_DIR}
  ${_mod_dir}/CheckStdVariant.cpp
  CXX_STANDARD ${FUB_CORE_CXX_STANDARD})
if (_has_std_variant)
  message(STATUS "Checking for std::variant... Success.")
else()
  message(STATUS "Checking for std::variant... Failure.")
endif()

set(FUB_CORE_USE_STD_VARIANT OFF)
set(FUB_CORE_USE_STD_EXPERIMENTAL_VARIANT OFF)
if (_has_std_variant)
  set(FUB_CORE_USE_STD_VARIANT ON)
else()
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/variant/.git)
    execute_process(COMMAND git submodule update --init --remote --depth=1 -- third-party/variant
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  set(FUB_CORE_USE_MPARK_VARIANT ON)
  add_library(mpark_variant INTERFACE)
  target_include_directories(mpark_variant INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/third-party/variant/include)
endif()

mark_as_advanced(FUB_CORE_USE_STD_VARIANT)
mark_as_advanced(FUB_CORE_USE_STD_EXPERIMENTAL_VARIANT)
mark_as_advanced(FUB_CORE_USE_MPARK_VARIANT)

unset(_mod_dir)
unset(_has_std_variant)
unset(_has_std_experimental_variant)
