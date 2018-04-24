# This script checks whether an implementation for std::invoke is available.

get_filename_component(_mod_dir ${CMAKE_CURRENT_LIST_FILE} PATH)

message(STATUS "Checking for std::invoke...")
try_compile(_has_std_invoke
  ${CMAKE_BINARY_DIR}
  ${_mod_dir}/CheckStdInvoke.cpp
  CXX_STANDARD ${FUB_CORE_CXX_STANDARD})
message(STATUS "Checking for std::invoke... done.")

if (_has_std_invoke)
  set(FUB_CORE_USE_STD_INVOKE ON)
else()
  set(FUB_CORE_USE_STD_INVOKE OFF)
endif()
mark_as_advanced(FUB_CORE_USE_STD_INVOKE)

unset(_mod_dir)
unset(_has_std_invoke)
