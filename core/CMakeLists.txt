cmake_minimum_required(VERSION 3.10)

project(ReactiveEulerSolver_core CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

###############################################################################
#                                         Set the default C++ standard to C++14

if (NOT DEFINED FUB_WITH_CXX_STANDARD)
  set(FUB_WITH_CXX_STANDARD 14)
endif()

###############################################################################
#                                             Find Boost or throw a fatal error

find_package(Boost REQUIRED)

###############################################################################
#                                      Download fmtlib and add the subdirectory

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/fmt/.git)
  execute_process(COMMAND git submodule update --remote --init --depth=1 -- third-party/fmt
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/third-party/fmt)


###############################################################################
#                                                Find Vc or throw a fatal error

find_package(Vc REQUIRED)

###############################################################################
#                                              Add range-v3 and define a target

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

configure_file(include/fub/config.hpp.in ${CMAKE_BINARY_DIR}/fub/core/config.hpp)

add_library(ReactiveEulerSolver_core INTERFACE)

# Add config.hpp to include directories
target_include_directories(ReactiveEulerSolver_core
  INTERFACE include ${CMAKE_BINARY_DIR} ${Vc_INCLUDE_DIR})

# Always Add dependencies to boost range-v3 and fmtlib
target_link_libraries(ReactiveEulerSolver_core
  INTERFACE Boost::boost range-v3 fmt Vc)

# If std::optional is not available we add a dependency to nonstd::optional
if (NOT FUB_WITH_STD_OPTIONAL)
  target_link_libraries(ReactiveEulerSolver_core INTERFACE optional-lite)
endif()

# If std::variant is not available we add a dependency to mpark::variant
if (NOT FUB_WITH_STD_VARIANT)
  target_link_libraries(ReactiveEulerSolver_core INTERFACE mpark_variant)
endif()

# Set the minmum required C++ to $FUB_WITH_CXX_STANDARD
set_target_properties(ReactiveEulerSolver_core
  PROPERTIES INTERFACE_COMPILE_FEATURES cxx_std_${FUB_WITH_CXX_STANDARD})

option(FUB_CORE_WITH_TESTS "Whether to compile tests." ON)

if (FUB_CORE_WITH_TESTS)
  add_subdirectory(tests)
endif()
