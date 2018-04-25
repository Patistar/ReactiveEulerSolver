# This script checks whether an implementation for std::apply is available.

get_filename_component(_mod_dir ${CMAKE_CURRENT_LIST_FILE} PATH)

message(STATUS "Checking for std::apply...")
try_compile(_has_std_apply
  ${CMAKE_BINARY_DIR}
  ${_mod_dir}/CheckStdApply.cpp
  CXX_STANDARD ${FUB_CORE_CXX_STANDARD})
if (_has_std_apply)
  message(STATUS "Checking for std::apply... Success")
else()
  message(STATUS "Checking for std::apply... Failure")
endif()

if (_has_std_apply)
  set(FUB_CORE_USE_STD_APPLY ON)
else()
  set(FUB_CORE_USE_STD_APPLY OFF)
endif()
mark_as_advanced(FUB_CORE_USE_STD_APPLY)

unset(_mod_dir)
unset(_has_std_apply)
