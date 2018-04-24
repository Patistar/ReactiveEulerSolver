cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_apps CXX)

if (NOT EXISTS core/.git)
  execute_process(COMMAND git submodule update --remote --init --depth=1 -- core
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(core)
endif()

if (NOT EXISTS grid/.git)
  execute_process(COMMAND git submodule update --remote --init --depth=1 -- grid
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(grid)
endif()

if (NOT EXISTS solver/.git)
  execute_process(COMMAND git submodule update --remote --init --depth=1 -- solver
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(solver)
endif()

if (NOT EXISTS driver/.git)
  execute_process(COMMAND git submodule update --remote --init --depth=1 -- driver
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory(driver)
endif()

add_subdirectory(single_stage)
add_subdirectory(burke_2012)
add_subdirectory(gri_30)
