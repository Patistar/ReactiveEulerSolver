# This script checks whether an implementation for std::optional is available.

include(CheckCXXSourceCompiles)

get_filename_component(_mod_dir ${CMAKE_CURRENT_LIST_FILE} PATH)

set(CMAKE_REQUIRED_FLAGS "-std=c++${FUB_CORE_CXX_STANDARD}")
check_cxx_source_compiles("
#include <optional>
int main() { std::optional<int> _; }
"
        FUB_CORE_USE_STD_OPTIONAL)


if (FUB_CORE_USE_STD_OPTIONAL)
  add_library(ReactiveEulerSolver_core_optional INTERFACE)
  target_compile_features(ReactiveEulerSolver_core_optional INTERFACE cxx_std_17)
else()
  set(CMAKE_REQUIRED_FLAGS "-std=c++${FUB_CORE_CXX_STANDARD}")
  check_cxx_source_compiles("
  #include <experimental/optional>
  int main() { std::experimental::optional<int> _; }
  "
          FUB_CORE_USE_STD_EXPERIMENTAL_OPTIONAL)
  if (FUB_CORE_USE_STD_EXPERIMENTAL_OPTIONAL)
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
