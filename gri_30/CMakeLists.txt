cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_apps_gri_30 CXX Fortran)

add_executable(serial.kinetic.gri_30.sod serial.kinetic.gri_30.sod.cpp)
target_link_libraries(serial.kinetic.gri_30.sod ReactiveEulerSolver.driver.serial.kinetic.gri_30.1d FortranRadauSolver)

find_package(HPX)

if (HPX_FOUND)
add_executable(hpx.kinetic.gri_30.sod hpx.kinetic.gri_30.sod.cpp)
target_link_libraries(hpx.kinetic.gri_30.sod ReactiveEulerSolver.driver.hpx.kinetic.gri_30.1d FortranRadauSolver)
hpx_setup_target(hpx.kinetic.gri_30.sod)

add_executable(hpx.kinetic.gri_30.detonation hpx.kinetic.gri_30.detonation.cpp)
target_link_libraries(hpx.kinetic.gri_30.detonation ReactiveEulerSolver.driver.hpx.kinetic.gri_30.1d FortranRadauSolver)
hpx_setup_target(hpx.kinetic.gri_30.detonation)
endif()
