cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_solver CXX Fortran)

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
