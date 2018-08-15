cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_solver CXX Fortran)

add_subdirectory(third-party/blaze)

add_library(FortranRadauSolver
  source/fub/ode_solver/radau.f
  source/fub/ode_solver/dc_decsol.f
  source/fub/ode_solver/decsol.f)

add_library(ReactiveEulerSolver.solver 
  source/fub/ode_solver/Radau.cpp)
target_include_directories(ReactiveEulerSolver.solver PUBLIC include PRIVATE source)
target_link_libraries(ReactiveEulerSolver.solver PUBLIC ReactiveEulerSolver.grid FortranRadauSolver blaze)

add_subdirectory(tests)
