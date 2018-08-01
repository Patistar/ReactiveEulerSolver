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
# FindP4EST
# --------
#
# Find the native P4EST include directories and libraries.
#
# Usage:
#
#   find_package(P4EST [major.[minor]] [EXACT] [QUIET] [REQUIRED])
#
# Users may modify the behaviour of this module with the following variables:
#
# * P4EST_DIR
# * P4EST_INCLUDE_DIR
# * P4EST_LIBRARY
#
# The module will set the following variables:
#
# * P4EST_FOUND
#
# This module will also create the "p4est" target thay may be used when building
# executables and libraries.

include(FindPackageHandleStandardArgs)

if (NOT P4EST_FOUND)
  set(P4EST_SEARCH_DIR ${P4EST_DIR})

  find_path(P4EST_INCLUDE_DIRS p4est_config.h
    HINTS ${P4EST_INCLUDE_DIR} ${P4EST_SEARCH_DIR}
    PATHS ${P4EST_DEFAULT_SEARCH_DIR}
    PATH_SUFFIXES include)

  find_library(P4EST_LIBRARIES p4est
    HINTS ${P4EST_LIBRARY} ${P4EST_SEARCH_DIR}
    PATHS ${P4EST_DEFAULT_SEARCH_DIR} ENV LIBRARY_PATH
    PATH_SUFFIXES lib)

  find_library(SC_LIBRARIES sc
    HINTS ${P4EST_LIBRARY} ${P4EST_SEARCH_DIR}
    PATHS ${P4EST_DEFAULT_SEARCH_DIR} ENV LIBRARY_PATH
    PATH_SUFFIXES lib)

  find_package_handle_standard_args(SC
    REQUIRED_VARS P4EST_INCLUDE_DIRS SC_LIBRARIES
    HANDLE_COMPONENTS)

  find_package_handle_standard_args(P4EST
    REQUIRED_VARS P4EST_INCLUDE_DIRS P4EST_LIBRARIES
    HANDLE_COMPONENTS)

  if(NOT CMAKE_VERSION VERSION_LESS 3.0 AND P4EST_FOUND AND SC_FOUND)
    add_library(sc SHARED IMPORTED)
    set_target_properties(sc PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES  ${P4EST_INCLUDE_DIRS}
      IMPORTED_LOCATION              ${SC_LIBRARIES})

    add_library(p4est SHARED IMPORTED)
    set_target_properties(p4est PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES  ${P4EST_INCLUDE_DIRS}
      IMPORTED_LOCATION              ${P4EST_LIBRARIES}
      IMPORTED_LINK_LIBRARIES sc)
  endif()

endif()

