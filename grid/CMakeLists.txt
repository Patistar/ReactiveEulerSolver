cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver_grid CXX)

if (NOT TARGET ReactiveEulerSolver_core)
  message(FATAL_ERROR "You can not use this module without ReactiveEulerSolver_core.")
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Check for CGNS
find_package(CGNS REQUIRED)
if (TARGET cgns)
  list(APPEND FUB_GRID_DEPENDENCIES cgns)
  set(FUB_WITH_CGNS TRUE)
endif()

# Check for P4EST
find_package(P4EST REQUIRED)
if (TARGET p4est AND TARGET sc)
  list(APPEND FUB_GRID_DEPENDENCIES z)
  list(APPEND FUB_GRID_DEPENDENCIES sc)
  list(APPEND FUB_GRID_DEPENDENCIES p4est)
  set(FUB_WITH_P4EST TRUE)
endif()

# Check for target "Boost::container"

if (TARGET Boost::container)
  list(APPEND FUB_GRID_DEPENDENCIES Boost::container)
  set(FUB_WITH_POLYMORPHIC_ALLOCATOR TRUE)
else()
  find_package(Boost REQUIRED container)
  list(APPEND FUB_GRID_DEPENDENCIES Boost::container)
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

add_library(ReactiveEulerSolver.grid
  source/fub/output/cgns.cpp
  source/fub/p4est/quadrant.cpp
  source/fub/p4est/connectivity.cpp
  source/fub/p4est/ghost_layer.cpp
  source/fub/p4est/tree.cpp
  source/fub/p4est/mesh.cpp
  source/fub/p4est/forest.cpp)
target_include_directories(ReactiveEulerSolver.grid PUBLIC include)
target_link_libraries(ReactiveEulerSolver.grid PUBLIC 
  ReactiveEulerSolver_core ${FUB_GRID_DEPENDENCIES} mpi)

add_subdirectory(tests)
# add_subdirectory(example)