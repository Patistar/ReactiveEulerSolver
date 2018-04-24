cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_grid CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(CGNS)

if (NOT TARGET ReactiveEulerSolver_core)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/core/.git)
    execute_process(COMMAND git submodule update --init --depth=1 -- third-party/core
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  add_subdirectory(third-party/core)
endif()

if (NOT TARGET Boost::container)
  find_package(Boost 1.65 REQUIRED COMPONENTS container)
endif()

if (NOT TARGET eigen)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/eigen-git-mirror/.git)
    execute_process(COMMAND git submodule update --init --depth=1 -- third-party/eigen-git-mirror
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  add_library(eigen INTERFACE)
  target_include_directories(eigen
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/third-party/eigen-git-mirror)
endif()

if (CGNS_FOUND)
add_library(ReactiveEulerSolver.grid source/fub/output/cgns.cpp)
target_include_directories(ReactiveEulerSolver.grid PUBLIC include)
target_link_libraries(ReactiveEulerSolver.grid PUBLIC ReactiveEulerSolver_core Boost::container eigen cgns)
else()
add_library(ReactiveEulerSolver.grid INTERFACE)
target_include_directories(ReactiveEulerSolver.grid INTERFACE include)
target_link_libraries(ReactiveEulerSolver.grid INTERFACE ReactiveEulerSolver_core Boost::container eigen)
endif()

add_subdirectory(tests)
