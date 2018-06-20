cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_apps_burke_2012 CXX Fortran)

add_executable(serial.burke_2012.sod serial.burke_2012.sod.cpp)
target_link_libraries(serial.burke_2012.sod 
  ReactiveEulerSolver.driver.serial.burke_2012.1d)

add_executable(serial.kinetic.burke_2012.sod serial.kinetic.burke_2012.sod.cpp)
target_link_libraries(serial.kinetic.burke_2012.sod
  ReactiveEulerSolver.driver.serial.kinetic.burke_2012.1d
  FortranRadauSolver)

add_executable(serial.kinetic.burke_2012.detonation serial.kinetic.burke_2012.detonation.cpp)
target_link_libraries(serial.kinetic.burke_2012.detonation
  ReactiveEulerSolver.driver.serial.kinetic.burke_2012.1d)

find_package(HPX)

if (HPX_FOUND)
add_executable(parallel.burke_2012.sod parallel.burke_2012.sod.cpp)
target_link_libraries(parallel.burke_2012.sod ReactiveEulerSolver.driver.parallel.burke_2012.1d)
hpx_setup_target(parallel.burke_2012.sod)

add_executable(parallel.kinetic.burke_2012.sod parallel.kinetic.burke_2012.sod.cpp)
target_link_libraries(parallel.kinetic.burke_2012.sod ReactiveEulerSolver.driver.parallel.kinetic.burke_2012.1d)
hpx_setup_target(parallel.kinetic.burke_2012.sod)

add_executable(parallel.kinetic.burke_2012.detonation parallel.kinetic.burke_2012.detonation.cpp)
target_link_libraries(parallel.kinetic.burke_2012.detonation ReactiveEulerSolver.driver.parallel.kinetic.burke_2012.1d)
hpx_setup_target(parallel.kinetic.burke_2012.detonation)
endif()

