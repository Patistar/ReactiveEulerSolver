cmake_minimum_required (VERSION 3.10)

project(ReactiveEulerSolver.grid_tests CXX)

add_executable(test.grid.variable_list test.grid.variable_list.cpp)
target_link_libraries(test.grid.variable_list ReactiveEulerSolver.grid Catch2)

add_executable(test.grid.variable_data test.grid.variable_data.cpp)
target_link_libraries(test.grid.variable_data ReactiveEulerSolver.grid Catch2)

add_executable(test.grid.variable_view test.grid.variable_view.cpp)
target_link_libraries(test.grid.variable_view ReactiveEulerSolver.grid Catch2)

add_executable(test.grid.gnuplot test.grid.gnuplot.cpp)
target_link_libraries(test.grid.gnuplot ReactiveEulerSolver.grid)

add_executable(test.grid.cgns test.grid.cgns.cpp)
target_link_libraries(test.grid.cgns ReactiveEulerSolver.grid)

add_executable(test.grid.p4est test.grid.p4est.cpp)
target_link_libraries(test.grid.p4est ReactiveEulerSolver.grid)

add_executable(test.grid.uniform_cartesian_coordinates 
  test.grid.uniform_cartesian_coordinates.cpp)
target_link_libraries(test.grid.uniform_cartesian_coordinates 
  ReactiveEulerSolver.grid Catch2)
