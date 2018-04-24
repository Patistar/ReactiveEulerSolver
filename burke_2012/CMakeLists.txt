cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_apps_burke_2012 CXX Fortran)

add_executable(serial.burke_2012.sod serial.burke_2012.sod.cpp)
target_link_libraries(serial.burke_2012.sod ReactiveEulerSolver.driver.serial.burke_2012.1d)

add_executable(serial.kinetic.burke_2012.sod serial.kinetic.burke_2012.sod.cpp)
target_link_libraries(serial.kinetic.burke_2012.sod ReactiveEulerSolver.driver.serial.kinetic.burke_2012.1d FortranRadauSolver)

add_executable(serial.kinetic.burke_2012.detonation serial.kinetic.burke_2012.detonation.cpp)
target_link_libraries(serial.kinetic.burke_2012.detonation ReactiveEulerSolver.driver.serial.kinetic.burke_2012.1d FortranRadauSolver)

find_package(HPX)

if (HPX_FOUND)
add_executable(hpx.burke_2012.sod hpx.burke_2012.sod.cpp)
target_link_libraries(hpx.burke_2012.sod ReactiveEulerSolver.driver.hpx.burke_2012.1d)
hpx_setup_target(hpx.burke_2012.sod)

add_executable(hpx.kinetic.burke_2012.sod hpx.kinetic.burke_2012.sod.cpp)
target_link_libraries(hpx.kinetic.burke_2012.sod ReactiveEulerSolver.driver.hpx.kinetic.burke_2012.1d FortranRadauSolver)
hpx_setup_target(hpx.kinetic.burke_2012.sod)

add_executable(hpx.kinetic.burke_2012.detonation hpx.kinetic.burke_2012.detonation.cpp)
target_link_libraries(hpx.kinetic.burke_2012.detonation ReactiveEulerSolver.driver.hpx.kinetic.burke_2012.1d FortranRadauSolver)
hpx_setup_target(hpx.kinetic.burke_2012.detonation)
endif()

