# ORIGINALLY FROM: https://github.com/Kitware/VTK/blob/master/CMake/FindTBB.cmake
# MODIFICATIONS BY orlandini
# - Find ThreadingBuildingBlocks include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(TBB
#    [REQUIRED]             # Fail with error if TBB is not found
#    )                      #
# Once done, this will define
#
#  TBB_FOUND - system has TBB
#  TBB_INCLUDE_DIRS - the TBB include directories
#  TBB_LIBRARIES - TBB libraries to be lined, doesn't include malloc or
#                  malloc proxy
#  TBB::tbb - imported target for the TBB library
#
#  TBB_VERSION_MAJOR - Major Product Version Number
#  TBB_VERSION_MINOR - Minor Product Version Number
#  TBB_INTERFACE_VERSION - Engineering Focused Version Number
#  TBB_COMPATIBLE_INTERFACE_VERSION - The oldest major interface version
#                                     still supported. This uses the engineering
#                                     focused interface version numbers.
#
#  TBB_MALLOC_FOUND - system has TBB malloc library
#  TBB_MALLOC_INCLUDE_DIRS - the TBB malloc include directories
#  TBB_MALLOC_LIBRARIES - The TBB malloc libraries to be lined
#  TBB::malloc - imported target for the TBB malloc library
#
#  TBB_MALLOC_PROXY_FOUND - system has TBB malloc proxy library
#  TBB_MALLOC_PROXY_INCLUDE_DIRS = the TBB malloc proxy include directories
#  TBB_MALLOC_PROXY_LIBRARIES - The TBB malloc proxy libraries to be lined
#  TBB::malloc_proxy - imported target for the TBB malloc proxy library
#
#
# This module reads hints about search locations from variables:
#  ENV TBB_ARCH_PLATFORM - for eg. set it to "mic" for Xeon Phi builds
#  ENV TBB_ROOT or just TBB_ROOT - root directory of tbb installation
#  ENV TBB_BUILD_PREFIX - specifies the build prefix for user built tbb
#                         libraries. Should be specified with ENV TBB_ROOT
#                         and optionally...
#  ENV TBB_BUILD_DIR - if build directory is different than ${TBB_ROOT}/build
#
#
# Modified by Robert Maynard from the original OGRE source
#
#-------------------------------------------------------------------
# This file is part of the CMake build system for OGRE
#     (Object-oriented Graphics Rendering Engine)
# For the latest info, see http://www.ogre3d.org/
#
# The contents of this file are placed in the public domain. Feel
# free to make use of it in any way you like.
#-------------------------------------------------------------------
#
# =========================================================================
# Taken from Copyright.txt in the root of the VTK source tree as per
# instructions to substitute the full license in place of the summary
# reference when distributing outside of VTK
# =========================================================================
#
#  Program:   Visualization Toolkit
#  Module:    Copyright.txt
#
# Copyright (c) 1993-2015 Ken Martin, Will Schroeder, Bill Lorensen
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither name of Ken Martin, Will Schroeder, or Bill Lorensen nor the names
#   of any contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# =========================================================================*/

#=============================================================================
#  FindTBB helper functions and macros
#

#====================================================
# Fix the library path in case it is a linker script
#====================================================
function(tbb_extract_real_library library real_library)
  if(NOT UNIX OR NOT EXISTS ${library})
    set(${real_library} "${library}" PARENT_SCOPE)
    return()
  endif()

  #Read in the first 4 bytes and see if they are the ELF magic number
  set(_elf_magic "7f454c46")
  file(READ ${library} _hex_data OFFSET 0 LIMIT 4 HEX)
  if(_hex_data STREQUAL _elf_magic)
    #we have opened a elf binary so this is what
    #we should link to
    set(${real_library} "${library}" PARENT_SCOPE)
    return()
  endif()

  file(READ ${library} _data OFFSET 0 LIMIT 1024)
  if("${_data}" MATCHES "INPUT \\(([^(]+)\\)")
    #extract out the .so name from REGEX MATCH command
    set(_proper_so_name "${CMAKE_MATCH_1}")

    #construct path to the real .so which is presumed to be in the same directory
    #as the input file
    get_filename_component(_so_dir "${library}" DIRECTORY)
    set(${real_library} "${_so_dir}/${_proper_so_name}" PARENT_SCOPE)
  else()
    #unable to determine what this library is so just hope everything works
    #and pass it unmodified.
    set(${real_library} "${library}" PARENT_SCOPE)
  endif()
endfunction()

#===============================================
# Do the final processing for the package find.
#===============================================
macro(findpkg_finish PREFIX TARGET_NAME)
  if (${PREFIX}_INCLUDE_DIR AND ${PREFIX}_LIBRARY)
    set(${PREFIX}_FOUND TRUE)
    set (${PREFIX}_INCLUDE_DIRS ${${PREFIX}_INCLUDE_DIR})
    set (${PREFIX}_LIBRARIES ${${PREFIX}_LIBRARY})
  else ()
    if (${PREFIX}_FIND_REQUIRED AND NOT ${PREFIX}_FIND_QUIETLY)
      message(FATAL_ERROR "Required library ${PREFIX} not found.")
    endif ()
  endif ()

  if (NOT TARGET "TBB::${TARGET_NAME}")
    if (${PREFIX}_LIBRARY_RELEASE)
      tbb_extract_real_library(${${PREFIX}_LIBRARY_RELEASE} real_release)
    endif ()
    if (${PREFIX}_LIBRARY_DEBUG)
      tbb_extract_real_library(${${PREFIX}_LIBRARY_DEBUG} real_debug)
    endif ()
    add_library(TBB::${TARGET_NAME} INTERFACE IMPORTED)
    set_target_properties(TBB::${TARGET_NAME} PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${${PREFIX}_INCLUDE_DIR}")
    if (${PREFIX}_LIBRARY_DEBUG AND ${PREFIX}_LIBRARY_RELEASE)
      set_target_properties(TBB::${TARGET_NAME} PROPERTIES
        INTERFACE_LINK_LIBRARIES "${real_release}"
        INTERFACE_LINK_LIBRARIES_DEBUG "${real_debug}"
        INTERFACE_LINK_LIBRARIES_RELEASE "${real_release}")
    elseif (${PREFIX}_LIBRARY_RELEASE)
      set_target_properties(TBB::${TARGET_NAME} PROPERTIES
        INTERFACE_LINK_LIBRARIES "${real_release}")
    elseif (${PREFIX}_LIBRARY_DEBUG)
      set_target_properties(TBB::${TARGET_NAME} PROPERTIES
        INTERFACE_LINK_LIBRARIES "${real_debug}")
    endif ()
  endif ()

  #mark the following variables as internal variables
  mark_as_advanced(${PREFIX}_INCLUDE_DIR
                   ${PREFIX}_LIBRARY
                   ${PREFIX}_LIBRARY_DEBUG
                   ${PREFIX}_LIBRARY_RELEASE)
endmacro()

#===============================================
# Generate debug names from given release names
#===============================================
macro(get_debug_names PREFIX)
  foreach(i ${${PREFIX}})
    set(${PREFIX}_DEBUG ${${PREFIX}_DEBUG} ${i}d ${i}D ${i}_d ${i}_D ${i}_debug ${i})
  endforeach()
endmacro()

#===============================================
# See if we have env vars to help us find tbb
#===============================================
macro(getenv_path VAR)
   set(ENV_${VAR} $ENV{${VAR}})
   # replace won't work if var is blank
   if (ENV_${VAR})
     string( REGEX REPLACE "\\\\" "/" ENV_${VAR} ${ENV_${VAR}} )
   endif ()
endmacro()

#===============================================
# Couple a set of release AND debug libraries
#===============================================
macro(make_library_set PREFIX)
    #@orlandini until i know better how to deal with is, picking release version if it is found
    if (NOT ${PREFIX}_RELEASE AND ${PREFIX}_DEBUG)
    set(${PREFIX} ${${PREFIX}_DEBUG})
  else()
    set(${PREFIX} ${${PREFIX}_RELEASE})  
  endif ()
endmacro()

#===============================================
# Ensure that the release & debug libraries found are from the same installation.
#===============================================
macro(find_tbb_library_verifying_release_debug_locations PREFIX)
  find_library(${PREFIX}_RELEASE
    NAMES ${_tbb_libname_prefix}${${PREFIX}_NAMES}${_tbb_shared_lib}
    HINTS ${_tbb_search_paths}
    PATH_SUFFIXES ${_tbb_libpath_suffix})
if(CMAKE_TBB_DEBUG)
    message(STATUS "${PREFIX} RELEASE ${${PREFIX}_RELEASE}")
endif(CMAKE_TBB_DEBUG)
  if (${PREFIX}_RELEASE)
    # To avoid finding a mismatched set of release & debug libraries from
    # different installations if the first found does not have debug libraries
    # by forcing the search for debug to only occur within the detected release
    # library directory (if found).  Although this would break detection if the
    # release & debug libraries were shipped in different directories, this is
    # not the case in the official TBB releases for any platform.
    get_filename_component(
      FOUND_RELEASE_LIB_DIR "${${PREFIX}_RELEASE}" DIRECTORY)
    find_library(${PREFIX}_DEBUG
      NAMES ${${PREFIX}_NAMES_DEBUG}
      HINTS ${FOUND_RELEASE_LIB_DIR}
      NO_DEFAULT_PATH)
  else()
    find_library(${PREFIX}_DEBUG
      NAMES ${${PREFIX}_NAMES_DEBUG}
      HINTS ${_tbb_search_paths}
      PATH_SUFFIXES ${_tbb_libpath_suffix})
  endif()
endmacro()

#=============================================================================
#  Now to actually find TBB
#

# Get path, convert backslashes as ${ENV_${var}}
getenv_path(TBB_ROOT)

# initialize search paths
set(TBB_PREFIX_PATH ${TBB_ROOT} ${ENV_TBB_ROOT})
set(_tbb_search_paths "")


# If user built from sources
set(TBB_BUILD_PREFIX $ENV{TBB_BUILD_PREFIX})
if (TBB_BUILD_PREFIX AND ENV_TBB_ROOT)
  getenv_path(TBB_BUILD_DIR)
  if (NOT ENV_TBB_BUILD_DIR)
    set(ENV_TBB_BUILD_DIR ${ENV_TBB_ROOT}/build)
  endif ()

  # include directory under ${ENV_TBB_ROOT}/include
  list(APPEND _tbb_search_paths
    ${ENV_TBB_BUILD_DIR}/${TBB_BUILD_PREFIX}_release
    ${ENV_TBB_BUILD_DIR}/${TBB_BUILD_PREFIX}_debug)
endif ()


# For Windows, let's assume that the user might be using the precompiled
# TBB packages from the main website. These use a rather awkward directory
# structure (at least for automatically finding the right files) depending
# on platform and compiler, but we'll do our best to accommodate it.
# Not adding the same effort for the precompiled linux builds, though. Those
# have different versions for CC compiler versions and linux kernels which
# will never adequately match the user's setup, so there is no feasible way
# to detect the "best" version to use. The user will have to manually
# select the right files. (Chances are the distributions are shipping their
# custom version of tbb, anyway, so the problem is probably nonexistent.)
if (WIN32 AND MSVC)
  set(COMPILER_PREFIX "vc7.1")
  if (MSVC_VERSION EQUAL 1400)
    set(COMPILER_PREFIX "vc8")
  elseif(MSVC_VERSION EQUAL 1500)
    set(COMPILER_PREFIX "vc9")
  elseif(MSVC_VERSION EQUAL 1600)
    set(COMPILER_PREFIX "vc10")
  elseif(MSVC_VERSION EQUAL 1700)
    set(COMPILER_PREFIX "vc11")
  elseif(MSVC_VERSION EQUAL 1800)
    set(COMPILER_PREFIX "vc12")
  elseif(MSVC_VERSION GREATER_EQUAL 1900)
    set(COMPILER_PREFIX "vc14")
  endif ()
endif ()

# For OS X binary distribution, choose libc++ based libraries for Mavericks (10.9)
# and above and AppleClang
if (CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND
    NOT CMAKE_SYSTEM_VERSION VERSION_LESS 13.0)
  set (USE_LIBCXX OFF)
  cmake_policy(GET CMP0025 POLICY_VAR)

  if (POLICY_VAR STREQUAL "NEW")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
      set (USE_LIBCXX ON)
    endif ()
  else ()
    if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      set (USE_LIBCXX ON)
    endif ()
  endif ()

  if (USE_LIBCXX)
    foreach (dir IN LISTS TBB_PREFIX_PATH)
      list (APPEND _tbb_search_paths ${dir}/lib/libc++ ${dir}/libc++/lib)
    endforeach ()
  endif ()
endif ()

# check compiler ABI
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(COMPILER_PREFIX)
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.8)
    list(APPEND COMPILER_PREFIX "gcc4.8")
  endif()
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
    list(APPEND COMPILER_PREFIX "gcc4.7")
  endif()
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.4)
    list(APPEND COMPILER_PREFIX "gcc4.4")
  endif()
  list(APPEND COMPILER_PREFIX "gcc4.1")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(COMPILER_PREFIX)
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.0) # Complete guess
    list(APPEND COMPILER_PREFIX "gcc4.8")
  endif()
  if (NOT CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.6)
    list(APPEND COMPILER_PREFIX "gcc4.7")
  endif()
  list(APPEND COMPILER_PREFIX "gcc4.4")
else() # Assume compatibility with 4.4 for other compilers
  list(APPEND COMPILER_PREFIX "gcc4.4")
endif ()

if (WIN32)
	set(PROGRAM_FILE_ENVVAR "PROGRAMFILES(x86)")
	file(TO_CMAKE_PATH "$ENV{${PROGRAM_FILE_ENVVAR}}" PRG_FOLD)
	list(APPEND _tbb_search_paths "${PRG_FOLD}/Intel/Composer XE/tbb") # default until ParallelStudioXE2015
	list(APPEND _tbb_search_paths "${PRG_FOLD}/IntelSWTools/compilers_and_libraries/windows/tbb") # default for ParallelStudioXE2016 and later
	list(APPEND _tbb_search_paths "${PRG_FOLD}/Intel/oneAPI/tbb") # default for oneAPI (2020 and later)
elseif(UNIX)#we need to test if we need sth different for apple
	foreach (_MKL_VER ${_MKL_TEST_VERSIONS})
		list(APPEND _tbb_search_paths "/opt/intel/composerxe-${_MKL_VER}/tbb") # default until ParallelStudioXE2015 (root permissions)
		list(APPEND _tbb_search_paths "$ENV{HOME}/intel/composerxe-${_MKL_VER}/tbb") # default until ParallelStudioXE2015 (no root permissions)
	endforeach()
    if(APPLE)
	    list(APPEND _tbb_search_paths "/opt/intel/compilers_and_libraries/mac/tbb") # default for ParallelStudioXE2016 and later (root permissions)
    else()
        list(APPEND _tbb_search_paths "/opt/intel/compilers_and_libraries/linux/tbb") # default for ParallelStudioXE2016 and later (root permissions)
    endif()
	list(APPEND _tbb_search_paths "/opt/intel/oneapi/tbb") # default for oneAPI (2020) and later (root permissions)
    if(APPLE)
	    list(APPEND _tbb_search_paths "$ENV{HOME}/intel/compilers_and_libraries/mac/tbb") # default for ParallelStudioXE2016 and later (no root permissions)
    else()
        list(APPEND _tbb_search_paths "$ENV{HOME}/intel/compilers_and_libraries/linux/tbb") # default for ParallelStudioXE2016 and later (no root permissions)
    endif()
    list(APPEND _tbb_search_paths "$ENV{HOME}/intel/oneapi/tbb") # default for oneAPI (2020) and later (no root permissions)
endif()

if(CMAKE_TBB_DEBUG)
    message(STATUS "Searching for TBB in:")
    foreach (dir ${_tbb_search_paths})
        message(${dir})
    endforeach ()         
endif(CMAKE_TBB_DEBUG)

find_path(TBB_INCLUDE_DIR tbb/tbb.h
          HINTS ${_tbb_search_paths}
          PATH_SUFFIXES
          "latest/include" #newer (2021) oneAPI version
          "include" #when installed by apt-get on debian
          )
if(CMAKE_TBB_DEBUG)
    message(STATUS "TBB INCLUDE DIR ${TBB_INCLUDE_DIR}")
endif(CMAKE_TBB_DEBUG)

#if we haven't found TBB no point on going any further
if (NOT TBB_INCLUDE_DIR)
  return()
endif ()


set(_tbb_libpath_suffix_orig "lib")
set(_tbb_libpath_suffix ${_tbb_libpath_suffix_orig})
list(APPEND _tbb_libpath_suffix "latest/${_tbb_libpath_suffix_orig}")
  
if(WIN32)
    list(APPEND _tbb_libpath_suffix "${_tbb_libpath_suffix_orig}_win")
    list(APPEND _tbb_libpath_suffix "latest/${_tbb_libpath_suffix_orig}_win")
    set(_tbb_libname_prefix "")
    set(_tbb_shared_lib "_dll.lib")
    set(_tbb_static_lib ".lib")
elseif(APPLE)
    list(APPEND _tbb_libpath_suffix "lib")
    list(APPEND _tbb_libpath_suffix "${_tbb_libpath_suffix_orig}_mac")
    list(APPEND _tbb_libpath_suffix "latest/${_tbb_libpath_suffix_orig}_mac")
    set(_tbb_libname_prefix "lib")
    set(_tbb_shared_lib ".dylib")
    set(_tbb_static_lib ".a")
else() # LINUX
    list(APPEND _tbb_libpath_suffix "${_tbb_libpath_suffix_orig}_lin")
    list(APPEND _tbb_libpath_suffix "latest/${_tbb_libpath_suffix_orig}_lin")
    set(_tbb_libname_prefix "lib")
    set(_tbb_shared_lib ".so")
    set(_tbb_static_lib ".a")
endif()

if(CMAKE_SIZEOF_VOID_P EQUAL 4) # 32 bit
    set(_tbb_arch_suffix "ia32")
else() # 64 bit
    set(_tbb_arch_suffix "intel64")
endif()


foreach (pre ${COMPILER_PREFIX})
  list(APPEND _tbb_libpath_suffix "latest/${_tbb_libpath_suffix_orig}/${pre}")
  list(APPEND _tbb_libpath_suffix "latest/${_tbb_libpath_suffix_orig}/${_tbb_arch_suffix}/${pre}")
endforeach ()

if(CMAKE_TBB_DEBUG)
    message(STATUS "Lib suffixes for TBB:")
    foreach (dir ${_tbb_libpath_suffix})
        message(${dir})
    endforeach ()         
endif(CMAKE_TBB_DEBUG)

set(TBB_LIBRARY_NAMES tbb)
get_debug_names(TBB_LIBRARY_NAMES)

find_tbb_library_verifying_release_debug_locations(TBB_LIBRARY)
make_library_set(TBB_LIBRARY)


findpkg_finish(TBB tbb)
#=============================================================================
# Look for TBB's malloc package
set(TBB_MALLOC_LIBRARY_NAMES tbbmalloc)
get_debug_names(TBB_MALLOC_LIBRARY_NAMES)

find_path(TBB_MALLOC_INCLUDE_DIR tbb/tbb.h
          HINTS ${_tbb_search_paths}
          PATH_SUFFIXES
          "latest/include" #newer (2021) oneAPI version
          "include" #when installed by apt-get on debian
          )
find_tbb_library_verifying_release_debug_locations(TBB_MALLOC_LIBRARY)
make_library_set(TBB_MALLOC_LIBRARY)

findpkg_finish(TBB_MALLOC tbbmalloc)

#=============================================================================
# Look for TBB's malloc proxy package
set(TBB_MALLOC_PROXY_LIBRARY_NAMES tbbmalloc_proxy)
get_debug_names(TBB_MALLOC_PROXY_LIBRARY_NAMES)

find_path(TBB_MALLOC_PROXY_INCLUDE_DIR tbb/tbbmalloc_proxy.h
          HINTS ${_tbb_search_paths}
          PATH_SUFFIXES
          "latest/include" #newer (2021) oneAPI version
          "include" #when installed by apt-get on debian
          )
find_tbb_library_verifying_release_debug_locations(TBB_MALLOC_PROXY_LIBRARY)
make_library_set(TBB_MALLOC_PROXY_LIBRARY)

findpkg_finish(TBB_MALLOC_PROXY tbbmalloc_proxy)


#=============================================================================
#parse all the version numbers from tbb
if(NOT TBB_VERSION)
  set(TBB_VERSION_FILE_PRIOR_TO_TBB_2021_1
    "${TBB_INCLUDE_DIR}/tbb/tbb_stddef.h")
  set(TBB_VERSION_FILE_AFTER_TBB_2021_1
    "${TBB_INCLUDE_DIR}/tbb/version.h")
  if (EXISTS "${TBB_VERSION_FILE_PRIOR_TO_TBB_2021_1}")
    set(TBB_VERSION_FILE "${TBB_VERSION_FILE_PRIOR_TO_TBB_2021_1}")
  elseif (EXISTS "${TBB_VERSION_FILE_AFTER_TBB_2021_1}")
    set(TBB_VERSION_FILE "${TBB_VERSION_FILE_AFTER_TBB_2021_1}")
  else()
    message(FATAL_ERROR "Found TBB installation: ${TBB_INCLUDE_DIR} "
      "missing version header.")
  endif()

 #only read the start of the file
 file(STRINGS
      "${TBB_VERSION_FILE}"
      TBB_VERSION_CONTENTS
      REGEX "VERSION")

  string(REGEX REPLACE
    ".*#define TBB_VERSION_MAJOR ([0-9]+).*" "\\1"
    TBB_VERSION_MAJOR "${TBB_VERSION_CONTENTS}")

  string(REGEX REPLACE
    ".*#define TBB_VERSION_MINOR ([0-9]+).*" "\\1"
    TBB_VERSION_MINOR "${TBB_VERSION_CONTENTS}")

  string(REGEX REPLACE
        ".*#define TBB_INTERFACE_VERSION ([0-9]+).*" "\\1"
        TBB_INTERFACE_VERSION "${TBB_VERSION_CONTENTS}")

  string(REGEX REPLACE
        ".*#define TBB_COMPATIBLE_INTERFACE_VERSION ([0-9]+).*" "\\1"
        TBB_COMPATIBLE_INTERFACE_VERSION "${TBB_VERSION_CONTENTS}")

endif()
