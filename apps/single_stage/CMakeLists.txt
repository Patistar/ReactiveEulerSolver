cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_apps_single_stage)

add_executable(serial.single_stage.sod serial.single_stage.sod.cpp)
target_link_libraries(serial.single_stage.sod ReactiveEulerSolver.driver.serial.single_stage.1d)

find_package(HPX)

if (HPX_FOUND)
add_executable(parallel.single_stage.sod parallel.single_stage.sod.cpp)
target_link_libraries(parallel.single_stage.sod ReactiveEulerSolver.driver.parallel.single_stage.1d)
hpx_setup_target(parallel.single_stage.sod)

add_executable(parallel.single_stage.2d.circle parallel.single_stage.2d.circle.cpp)
target_link_libraries(parallel.single_stage.2d.circle ReactiveEulerSolver.driver.parallel.single_stage.2d)
hpx_setup_target(parallel.single_stage.2d.circle)
endif()

