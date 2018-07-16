cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_grid CXX)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(CGNS REQUIRED)
list(APPEND FUB_GRID_DEPENDENCIES cgns)

if (NOT TARGET Boost::container)
  find_package(Boost COMPONENTS container)
  if (TARGET Boost::container)
    set(FUB_WITH_POLYMORPHIC_ALLOCATOR TRUE)
    list(APPEND FUB_GRID_DEPENDENCIES Boost::container)
  else()
    set(FUB_WITH_POLYMORPHIC_ALLOCATOR FALSE)
  endif()
else()
  set(FUB_WITH_POLYMORPHIC_ALLOCATOR TRUE)
endif()

configure_file(include/fub/config.hpp.in ${CMAKE_BINARY_DIR}/fub/grid/config.hpp)

if (NOT TARGET eigen)
  if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third-party/eigen-git-mirror/.git)
    execute_process(COMMAND git submodule update --init --depth=1 -- third-party/eigen-git-mirror
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
  add_library(eigen INTERFACE)
  target_include_directories(eigen
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/third-party/eigen-git-mirror)
endif()
list(APPEND FUB_GRID_DEPENDENCIES eigen)

add_library(ReactiveEulerSolver.grid_includes INTERFACE)
target_include_directories(ReactiveEulerSolver.grid_includes INTERFACE include)
target_link_libraries(ReactiveEulerSolver.grid_includes INTERFACE 
  ReactiveEulerSolver_core ${FUB_GRID_DEPENDENCIES})

add_library(ReactiveEulerSolver.grid source/fub/output/cgns.cpp)
target_include_directories(ReactiveEulerSolver.grid PUBLIC include)
target_link_libraries(ReactiveEulerSolver.grid PUBLIC 
  ReactiveEulerSolver_core ${FUB_GRID_DEPENDENCIES})

add_subdirectory(tests)
