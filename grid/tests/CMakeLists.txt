cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver.grid_tests CXX)

if (NOT TARGET Catch2)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/Catch2/.git)
    execute_process(COMMAND git submodule update --init -- third-party/Catch2
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  add_library(Catch2 INTERFACE)
  target_include_directories(Catch2
    INTERFACE
      ${CMAKE_CURRENT_SOURCE_DIR}/../third-party/Catch2/single_include)
endif()

add_executable(test.grid.patch test.grid.patch.cpp)
target_link_libraries(test.grid.patch ReactiveEulerSolver.grid Catch2)

add_executable(test.grid.octree test.grid.octree.cpp)
target_link_libraries(test.grid.octree ReactiveEulerSolver.grid Catch2)

if (NOT DEFINED HPX_FOUND OR NOT HPX_FOUND)
  find_package(HPX 1.0)
endif()

if (DEFINED HPX_FOUND AND HPX_FOUND)
  add_hpx_executable(test.grid.distributed_grid
    SOURCES test.grid.distributed_grid.cpp
    DEPENDENCIES ReactiveEulerSolver.grid Catch2)
endif()

add_executable(test.grid.variable_accessor test.grid.variable_accessor.cpp)
target_link_libraries(test.grid.variable_accessor ReactiveEulerSolver.grid_includes Catch2)

add_executable(test.grid.interval_map test.grid.interval_map.cpp)
target_link_libraries(test.grid.interval_map ReactiveEulerSolver.grid)

add_executable(test.grid.grid test.grid.grid.cpp)
target_link_libraries(test.grid.grid ReactiveEulerSolver.grid)

add_executable(test.grid.patch_view test.grid.patch_view.cpp)
target_link_libraries(test.grid.patch_view ReactiveEulerSolver.grid Catch2)

add_executable(test.grid.uniform_cartesian_coordinates test.grid.uniform_cartesian_coordinates.cpp)
target_link_libraries(test.grid.uniform_cartesian_coordinates ReactiveEulerSolver.grid Catch2)

add_executable(test.grid.variables test.grid.variables.cpp)
target_link_libraries(test.grid.variables ReactiveEulerSolver.grid Catch2)
