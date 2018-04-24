# This script checks whether an implementation for std::clamp is available.

get_filename_component(_mod_dir ${CMAKE_CURRENT_LIST_FILE} PATH)

message(STATUS "Checking for std::clamp...")
try_compile(_has_std_clamp
  ${CMAKE_BINARY_DIR}
  ${_mod_dir}/CheckStdClamp.cpp
  CXX_STANDARD ${FUB_CORE_CXX_STANDARD})
message(STATUS "Checking for std::clamp... done")

if (_has_std_clamp)
  set(FUB_CORE_USE_STD_CLAMP ON)
else()
  set(FUB_CORE_USE_STD_CLAMP OFF)
endif()
mark_as_advanced(FUB_CORE_USE_STD_CLAMP)

unset(_mod_dir)
unset(_has_std_clamp)
