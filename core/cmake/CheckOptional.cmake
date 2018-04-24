# This script checks whether an implementation for std::optional is available.

get_filename_component(_mod_dir ${CMAKE_CURRENT_LIST_FILE} PATH)
message(STATUS "Checking for std::optional...")
try_compile(_has_std_optional
  ${CMAKE_BINARY_DIR}
  ${_mod_dir}/CheckStdOptional.cpp
  CXX_STANDARD ${FUB_CORE_CXX_STANDARD})
message(STATUS "Checking for std::optional... done.")

set(FUB_CORE_USE_STD_OPTIONAL OFF)
set(FUB_CORE_USE_STD_EXPERIMENTAL_OPTIONAL OFF)
set(FUB_CORE_USE_BOOST_OPTIONAL OFF)
if (_has_std_optional)
  set(FUB_CORE_USE_STD_OPTIONAL ON)
  add_library(ReactiveEulerSolver_core_optional INTERFACE)
  target_compile_features(ReactiveEulerSolver_core_optional INTERFACE cxx_std_17)
else()
  message(STATUS "Checking for std::experimental::optional...")
  try_compile(_has_std_experimental_optional
    ${CMAKE_BINARY_DIR}
    ${_mod_dir}/CheckStdExperimentalOptional.cpp
    CXX_STANDARD ${FUB_CORE_CXX_STANDARD})
  message(STATUS "Checking for std::experimental::optional... done.")
  if (_has_std_experimental_optional)
    set(FUB_CORE_USE_STD_EXPERIMENTAL_OPTIONAL ON)
    add_library(ReactiveEulerSolver_core_optional INTERFACE)
    target_compile_features(ReactiveEulerSolver_core_optional INTERFACE cxx_std_14)
  elseif(NOT TARGET Boost::boost)
    find_package(Boost REQUIRED)
    set(FUB_CORE_USE_BOOST_OPTIONAL ON)
    add_library(ReactiveEulerSolver_core_optional INTERFACE)
    target_compile_features(ReactiveEulerSolver_core_optional INTERFACE cxx_std_11)
  endif()
endif()

unset(_mod_dir)
unset(_has_std_optional)
unset(_has_std_experimental_optional)
