cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_solver_tests CXX)

add_executable(test.solver.burke_2012 test.solver.burke_2012.cpp)
target_link_libraries(test.solver.burke_2012 ReactiveEulerSolver.solver FortranRadauSolver)

add_executable(test.solver.ode_solver test.solver.ode_solver.cpp)
target_link_libraries(test.solver.ode_solver ReactiveEulerSolver.solver FortranRadauSolver blaze)