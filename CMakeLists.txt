cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver CXX Fortran)

add_subdirectory(core)
add_subdirectory(grid)
add_subdirectory(solver)
add_subdirectory(driver)
add_subdirectory(apps)

execute_process(COMMAND git submodule foreach git reset --hard
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
execute_process(COMMAND git submodule foreach git pull origin master
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
