cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_driver CXX Fortran)

if (NOT TARGET ReactiveEulerSolver_core)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/core/.git)
    execute_process(COMMAND git submodule update --init -- third-party/core
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  add_subdirectory(third-party/core)
endif()

if (NOT TARGET ReactiveEulerSolver.grid)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/grid/.git)
    execute_process(COMMAND git submodule update --init -- third-party/grid
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  add_subdirectory(third-party/grid)
endif()

if (NOT TARGET ReactiveEulerSolver.solver)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/solver/.git)
    execute_process(COMMAND git submodule update --init -- third-party/solver
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  add_subdirectory(third-party/solver)
endif()

if (NOT TARGET Boost::program_options)
  find_package(Boost 1.65 REQUIRED COMPONENTS program_options)
endif()

if (NOT DEFINED HPX_FOUND OR NOT HPX_FOUND)
  find_package(HPX 1.0)
endif()

add_library(ReactiveEulerSolver.driver.serial.single_stage.1d source/fub/serial/single_stage_1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.serial.single_stage.1d PUBLIC include)
target_link_libraries(ReactiveEulerSolver.driver.serial.single_stage.1d
  PUBLIC ReactiveEulerSolver.solver Boost::program_options)

if (HPX_FOUND)
add_library(ReactiveEulerSolver.driver.hpx.single_stage.1d source/fub/hpx/single_stage_1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.hpx.single_stage.1d PUBLIC include)
target_link_libraries(ReactiveEulerSolver.driver.hpx.single_stage.1d ReactiveEulerSolver.solver Boost::program_options)
hpx_setup_target(ReactiveEulerSolver.driver.hpx.single_stage.1d)
endif()

add_library(ReactiveEulerSolver.driver.serial.burke_2012.1d source/fub/serial/burke_2012_1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.serial.burke_2012.1d PUBLIC include)
target_link_libraries(ReactiveEulerSolver.driver.serial.burke_2012.1d PUBLIC ReactiveEulerSolver.solver Boost::program_options)

add_library(ReactiveEulerSolver.driver.serial.kinetic.burke_2012.1d source/fub/serial/kinetic.burke_2012_1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.serial.kinetic.burke_2012.1d PUBLIC include)
target_link_libraries(ReactiveEulerSolver.driver.serial.kinetic.burke_2012.1d PUBLIC ReactiveEulerSolver.solver Boost::program_options)

add_library(ReactiveEulerSolver.driver.serial.kinetic.gri_30.1d source/fub/serial/kinetic.gri_30.1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.serial.kinetic.gri_30.1d PUBLIC include)
target_link_libraries(ReactiveEulerSolver.driver.serial.kinetic.gri_30.1d PUBLIC ReactiveEulerSolver.solver Boost::program_options)

if (HPX_FOUND)
add_library(ReactiveEulerSolver.driver.hpx.burke_2012.1d source/fub/hpx/burke_2012_1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.hpx.burke_2012.1d PUBLIC include)
target_link_libraries(ReactiveEulerSolver.driver.hpx.burke_2012.1d ReactiveEulerSolver.solver Boost::program_options)
hpx_setup_target(ReactiveEulerSolver.driver.hpx.burke_2012.1d)

add_library(ReactiveEulerSolver.driver.hpx.kinetic.burke_2012.1d source/fub/hpx/kinetic.burke_2012_1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.hpx.kinetic.burke_2012.1d PUBLIC include)
target_link_libraries(ReactiveEulerSolver.driver.hpx.kinetic.burke_2012.1d ReactiveEulerSolver.solver Boost::program_options)
hpx_setup_target(ReactiveEulerSolver.driver.hpx.kinetic.burke_2012.1d)

add_library(ReactiveEulerSolver.driver.hpx.kinetic.gri_30.1d source/fub/hpx/kinetic.gri_30.1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.hpx.kinetic.gri_30.1d PUBLIC include)
target_link_libraries(ReactiveEulerSolver.driver.hpx.kinetic.gri_30.1d ReactiveEulerSolver.solver Boost::program_options)
hpx_setup_target(ReactiveEulerSolver.driver.hpx.kinetic.gri_30.1d)
endif()
