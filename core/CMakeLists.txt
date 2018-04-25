cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_core CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if (NOT DEFINED FUB_CORE_CXX_STANDARD)
  set(FUB_CORE_CXX_STANDARD 17)
endif()

find_package(Boost REQUIRED)

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/fmt/.git)
  execute_process(COMMAND git submodule update --remote --init --depth=1 -- third-party/fmt
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/third-party/fmt)
endif()

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/range-v3/.git)
  execute_process(COMMAND git submodule update --remote --init --depth=1 -- third-party/range-v3
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
add_library(range-v3 INTERFACE)
target_include_directories(range-v3
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/third-party/range-v3/include)


include(CheckOptional)
include(CheckVariant)
include(CheckClamp)
include(CheckInvoke)
include(CheckApply)
include(CheckArray)

configure_file(include/fub/config.hpp.in ${CMAKE_BINARY_DIR}/fub/core/config.hpp)

add_library(ReactiveEulerSolver_core INTERFACE)
target_include_directories(ReactiveEulerSolver_core
  INTERFACE include ${CMAKE_BINARY_DIR})
target_link_libraries(ReactiveEulerSolver_core 
  INTERFACE Boost::boost range-v3 fmt)
if (DEFINED FUB_CORE_USE_MPARK_VARIANT)
  target_link_libraries(ReactiveEulerSolver_core INTERFACE mpark_variant)
endif()
set_target_properties(ReactiveEulerSolver_core
  PROPERTIES INTERFACE_COMPILE_FEATURES cxx_std_${FUB_CORE_CXX_STANDARD})

option(FUB_CORE_WITH_TESTS "Whether to compile tests." ON)

if (FUB_CORE_WITH_TESTS)
  add_subdirectory(tests)
endif()
