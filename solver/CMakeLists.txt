cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_solver CXX Fortran)

if (NOT TARGET ReactiveEulerSolver_core)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/core/.git)
    execute_process(COMMAND git submodule update --init -- third-party/core
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  add_subdirectory(third-party/core)
endif()

if (NOT TARGET ReactiveEulerSolver.grid)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/grid/.git)
    execute_process(COMMAND git submodule update --init -- third-party/grid
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  add_subdirectory(third-party/grid)
endif()

add_library(FortranRadauSolver
  source/fub/ode_solver/radau.f
  source/fub/ode_solver/dc_decsol.f
  source/fub/ode_solver/decsol.f)

add_library(ReactiveEulerSolver.solver
  source/fub/ode_solver/Radau.cpp
  source/fub/euler/hll_riemann_solver.cpp
  source/fub/euler/hlle_signal_velocity.cpp
  source/fub/time_integrator/forward_euler.cpp
  source/fub/euler/mechanism/burke_2012.cpp
  source/fub/euler/mechanism/gri_30.cpp)
target_include_directories(ReactiveEulerSolver.solver PUBLIC include)
target_link_libraries(ReactiveEulerSolver.solver PUBLIC ReactiveEulerSolver.grid)

add_subdirectory(tests)
