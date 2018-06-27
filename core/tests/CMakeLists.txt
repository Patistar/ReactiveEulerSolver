cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_core_tests CXX)

if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../third-party/Catch2/.git)
  execute_process(COMMAND git submodule update --remote --init --depth=1 -- third-party/Catch2
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../)
endif()
add_library(Catch2 INTERFACE)
target_include_directories(Catch2
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../third-party/Catch2/single_include)

add_executable(test.core.fancy_pointer test.core.fancy_poiner.cpp)
target_link_libraries(test.core.fancy_pointer ReactiveEulerSolver_core Catch2)

add_executable(test.core.extents test.core.extents.cpp)
target_link_libraries(test.core.extents ReactiveEulerSolver_core)

add_executable(test.core.simd test.core.simd.cpp)
target_link_libraries(test.core.simd ReactiveEulerSolver_core)

add_executable(test.core.optional test.core.optional.cpp)
target_link_libraries(test.core.optional ReactiveEulerSolver_core)

add_executable(test.core.variant test.core.variant.cpp)
target_link_libraries(test.core.variant ReactiveEulerSolver_core)

add_executable(test.core.span test.core.span.cpp)
target_link_libraries(test.core.span ReactiveEulerSolver_core Catch2)

add_executable(test.core.mdspan test.core.mdspan.cpp)
target_link_libraries(test.core.mdspan ReactiveEulerSolver_core Catch2)
