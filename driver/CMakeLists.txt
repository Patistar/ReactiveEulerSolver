cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_driver CXX Fortran)

if (NOT TARGET Boost::program_options)
  find_package(Boost REQUIRED COMPONENTS program_options)
endif()

if (NOT DEFINED HPX_FOUND OR NOT HPX_FOUND)
  find_package(HPX 1.0)
endif()

add_library(ReactiveEulerSolver.driver.serial.single_stage.1d source/fub/serial/single_stage.1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.serial.single_stage.1d PUBLIC include PRIVATE source)
target_link_libraries(ReactiveEulerSolver.driver.serial.single_stage.1d
  PUBLIC ReactiveEulerSolver.solver Boost::program_options)

if (HPX_FOUND)
add_library(ReactiveEulerSolver.driver.parallel.single_stage.1d source/fub/parallel/single_stage.1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.parallel.single_stage.1d PUBLIC include PRIVATE source)
target_link_libraries(ReactiveEulerSolver.driver.parallel.single_stage.1d ReactiveEulerSolver.solver Boost::program_options)
hpx_setup_target(ReactiveEulerSolver.driver.parallel.single_stage.1d)

add_library(ReactiveEulerSolver.driver.distributed.single_stage.1d source/fub/distributed/single_stage.1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.distributed.single_stage.1d PUBLIC include PRIVATE source)
target_link_libraries(ReactiveEulerSolver.driver.distributed.single_stage.1d ReactiveEulerSolver.solver Boost::program_options)
hpx_setup_target(ReactiveEulerSolver.driver.distributed.single_stage.1d)

add_library(ReactiveEulerSolver.driver.parallel.single_stage.2d source/fub/parallel/single_stage.2d.cpp)
target_include_directories(ReactiveEulerSolver.driver.parallel.single_stage.2d PUBLIC include PRIVATE source)
target_link_libraries(ReactiveEulerSolver.driver.parallel.single_stage.2d ReactiveEulerSolver.solver Boost::program_options)
hpx_setup_target(ReactiveEulerSolver.driver.parallel.single_stage.2d)
endif()

add_library(ReactiveEulerSolver.driver.serial.burke_2012.1d source/fub/serial/burke_2012.1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.serial.burke_2012.1d PUBLIC include PRIVATE source)
target_link_libraries(ReactiveEulerSolver.driver.serial.burke_2012.1d PUBLIC ReactiveEulerSolver.solver Boost::program_options)

add_library(ReactiveEulerSolver.driver.serial.kinetic.burke_2012.1d source/fub/serial/kinetic.burke_2012.1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.serial.kinetic.burke_2012.1d PUBLIC include PRIVATE source)
target_link_libraries(ReactiveEulerSolver.driver.serial.kinetic.burke_2012.1d PUBLIC ReactiveEulerSolver.solver Boost::program_options)

add_library(ReactiveEulerSolver.driver.serial.kinetic.gri_30.1d source/fub/serial/kinetic.gri_30.1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.serial.kinetic.gri_30.1d PUBLIC include PRIVATE source)
target_link_libraries(ReactiveEulerSolver.driver.serial.kinetic.gri_30.1d PUBLIC ReactiveEulerSolver.solver Boost::program_options)

if (HPX_FOUND)
add_library(ReactiveEulerSolver.driver.parallel.burke_2012.1d source/fub/parallel/burke_2012.1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.parallel.burke_2012.1d PUBLIC include PRIVATE source)
target_link_libraries(ReactiveEulerSolver.driver.parallel.burke_2012.1d ReactiveEulerSolver.solver Boost::program_options)
hpx_setup_target(ReactiveEulerSolver.driver.parallel.burke_2012.1d)

add_library(ReactiveEulerSolver.driver.parallel.kinetic.burke_2012.1d source/fub/parallel/kinetic.burke_2012.1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.parallel.kinetic.burke_2012.1d PUBLIC include PRIVATE source)
target_link_libraries(ReactiveEulerSolver.driver.parallel.kinetic.burke_2012.1d ReactiveEulerSolver.solver Boost::program_options)
hpx_setup_target(ReactiveEulerSolver.driver.parallel.kinetic.burke_2012.1d)

add_library(ReactiveEulerSolver.driver.parallel.kinetic.gri_30.1d source/fub/parallel/kinetic.gri_30.1d.cpp)
target_include_directories(ReactiveEulerSolver.driver.parallel.kinetic.gri_30.1d PUBLIC include PRIVATE source)
target_link_libraries(ReactiveEulerSolver.driver.parallel.kinetic.gri_30.1d ReactiveEulerSolver.solver Boost::program_options)
hpx_setup_target(ReactiveEulerSolver.driver.parallel.kinetic.gri_30.1d)
endif()
