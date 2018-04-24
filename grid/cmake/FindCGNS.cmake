# Copyright (c) 2018 Maikel Nadolski
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#
# FindCGNS
# --------
#
# Find the native CGNS include directories and libraries.
#
# Usage:
#
#   find_package(CGNS [major.[minor]] [EXACT] [QUIET] [REQUIRED])
#
# Users may modify the behaviour of this module with the following variables:
#
# * CGNS_ROOT_DIR
# * CGNS_INCLUDE_DIR
# * CGNS_LIBRARY
#
# The module will set the following variables:
#
# * CGNS_FOUND
#
# This module will also create the "cgns" target thay may be used when building
# executables and libraries.

include(FindPackageHandleStandardArgs)

if (NOT CGNS_FOUND)

  set(CGNS_SEARCH_DIR ${CGNS_ROOT_DIR})

  find_path(CGNS_INCLUDE_DIRS cgnslib.h
    HINTS ${CGNS_INCLUDE_DIR} ${CGNS_SEARCH_DIR}
    PATHS ${CGNS_DEFAULT_SEARCH_DIR}
    PATH_SUFFIXES include)

  find_library(CGNS_LIBRARIES cgns
    HINTS ${CGNS_LIBRARY} ${CGNS_SEARCH_DIR}
    PATHS ${CGNS_DEFAULT_SEARCH_DIR} ENV LIBRARY_PATH
    PATH_SUFFIXES lib)

  find_package_handle_standard_args(CGNS
    REQUIRED_VARS CGNS_INCLUDE_DIRS CGNS_LIBRARIES
    HANDLE_COMPONENTS
    VERSION_VAR TBB_VERSION)

  if(NOT CMAKE_VERSION VERSION_LESS 3.0 AND CGNS_FOUND)
    add_library(cgns SHARED IMPORTED)
    set_target_properties(cgns PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES  ${CGNS_INCLUDE_DIRS}
      IMPORTED_LOCATION              ${CGNS_LIBRARIES})
  endif()

endif()

