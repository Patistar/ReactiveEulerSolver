get_filename_component(_mod_dir ${CMAKE_CURRENT_LIST_FILE} PATH)
message(STATUS "Checking for constexpr-enabled std::array...")
try_compile(_has_constexpr_std_array
  ${CMAKE_BINARY_DIR}
  ${_mod_dir}/CheckStdArray.cpp
  CXX_STANDARD ${FUB_CORE_CXX_STANDARD})
message(STATUS "Checking for constexpr-enabled std::array... done.")

set(FUB_CORE_USE_STD_ARRAY _has_constexpr_std_array)
mark_as_advanced(FUB_CORE_USE_STD_ARRAY)

unset(_mod_dir)
unset(_has_std_variant)
unset(_has_std_experimental_variant)
