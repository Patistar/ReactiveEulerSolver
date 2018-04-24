cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_apps CXX)

add_subdirectory(core)
add_subdirectory(grid)
add_subdirectory(solver)
add_subdirectory(driver)

add_subdirectory(single_stage)
add_subdirectory(burke_2012)
add_subdirectory(gri_30)
