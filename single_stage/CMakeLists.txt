cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_apps_single_stage)

add_executable(serial.single_stage.sod serial.single_stage.sod.cpp)
target_link_libraries(serial.single_stage.sod ReactiveEulerSolver.driver.serial.single_stage.1d)

find_package(HPX)

if (HPX_FOUND)
add_executable(hpx.single_stage.sod hpx.single_stage.sod.cpp)
target_link_libraries(hpx.single_stage.sod ReactiveEulerSolver.driver.hpx.single_stage.1d)
hpx_setup_target(hpx.single_stage.sod)
endif()

