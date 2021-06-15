#[[
The following CMake module has been based in two different FindMKL.cmake modules.

First, it was based on
https://gitlab.kitware.com/cmake/cmake/-/blob/34b308ec0b9616cdf273d697f0ef78f335c28511/Modules/FindMKL.cmake
by Teodor Nikolov

and subsequentely the search paths idea has been taken from:
https://github.com/projectchrono/chrono/blob/f33021edaad37746c3093150c8e3c69e1f000074/cmake/FindMKL.cmake
by Radu Serban and Dario Mongoni
]]
#[=======================================================================[.rst:
FindMKL
-------

The following conventions are used:

intel / INTEL  - Bindings for everything except GNU Fortran
gf / GF        - GNU Fortran bindings
seq / SEQ      - sequential MKL
omp / OMP      - threaded MKL with OpenMP back end
tbb / TBB      - threaded MKL with TBB back end
32bit / 32BIT  - MKL 32 bit integer interface (used most often)
64bit / 64BIT  - MKL 64 bit integer interface
mpich / MPICH  - MPICH / IntelMPI BLACS back end
ompi / OMPI    - OpenMPI BLACS back end
st / ST        - static libraries
dyn / DYN      - dynamic libraries

The module attempts to define a target for each MKL configuration. The
configuration will not be available if there are missing library files or a
missing dependency.

MKL Link line advisor:
  https://software.intel.com/content/www/us/en/develop/documentation/get-started-with-intel-oneapi-base-linux/top/before-you-begin.html
  https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

Note: Mixing GCC and Intel OpenMP backends is a bad idea.

Search variables
^^^^^^^^^^^^^^^^

``MKLROOT``
  Environment variable set to MKL's root directory

``MKL_ROOT``
  CMake variable set to MKL's root directory

Example usage
^^^^^^^^^^^^^

To Find MKL:

  find_package(MKL REQUIRED)

To check if target is available:

  if (TARGET mkl::scalapack_mpich_intel_32bit_omp_dyn)
    ...
  endif()

To link to an available target (see list below):

  target_link_libraries(... mkl::scalapack_mpich_intel_32bit_omp_dyn)

Note: dependencies are handled for you (MPI, OpenMP, ...)

Imported targets
^^^^^^^^^^^^^^^^

MKL (BLAS, LAPACK, FFT) targets:

  mkl::mkl_[gf|intel]_[32bit|64bit]_[seq|omp|tbb]_[st|dyn] e.g.

  mkl::mkl_intel_32bit_omp_dyn

BLACS targets:

  mkl::blacs_[mpich|ompi]_[gf|intel]_[32bit|64bit]_[seq|omp|tbb]_[st|dyn] e.g.

  mkl::blacs_intel_mpich_32bit_seq_st

ScaLAPACK targets:

  mkl::scalapack_[mpich|ompi]_[gf|intel]_[32bit|64bit]_[seq|omp|tbb]_[st|dyn] e.g.

  mkl::scalapack_intel_mpich_64bit_omp_dyn

Result variables
^^^^^^^^^^^^^^^^

MKL_FOUND

Not supported
^^^^^^^^^^^^^

- F95 interfaces

#]=======================================================================]

cmake_minimum_required(VERSION 3.12)

# Modules
#
include(FindPackageHandleStandardArgs)

if(NOT (CMAKE_C_COMPILER_LOADED OR
        CMAKE_CXX_COMPILER_LOADED OR
        CMAKE_Fortran_COMPILER_LOADED))
    message(FATAL_ERROR "FindMKL requires Fortran, C, or C++ to be enabled.")
endif()

# Dependencies
find_package(TBB QUIET)
find_package(Threads QUIET)
find_package(MPI COMPONENTS CXX QUIET)
find_package(OpenMP COMPONENTS CXX QUIET)

# If MKL_ROOT is not set, set it via the env variable MKLROOT.
#
if(NOT DEFINED MKL_ROOT)
    set(MKL_ROOT $ENV{MKLROOT} CACHE PATH "MKL's root directory.")
endif()

set(_mkl_search_paths ${MKL_ROOT})

# Add the default install location to the search path
if (WIN32)
	set(PROGRAM_FILE_ENVVAR "PROGRAMFILES(x86)")
	file(TO_CMAKE_PATH "$ENV{${PROGRAM_FILE_ENVVAR}}" PRG_FOLD)
	list(APPEND _mkl_search_paths "${PRG_FOLD}/Intel/Composer XE/mkl") # default until ParallelStudioXE2015
	list(APPEND _mkl_search_paths "${PRG_FOLD}/IntelSWTools/compilers_and_libraries/windows/mkl") # default for ParallelStudioXE2016 and later
	list(APPEND _mkl_search_paths "${PRG_FOLD}/Intel/oneAPI/mkl") # default for oneAPI (2020 and later)
elseif(UNIX)#we need to test if we need sth different for apple
	foreach (_MKL_VER ${_MKL_TEST_VERSIONS})
		list(APPEND _mkl_search_paths "/opt/intel/composerxe-${_MKL_VER}/mkl") # default until ParallelStudioXE2015 (root permissions)
		list(APPEND _mkl_search_paths "$ENV{HOME}/intel/composerxe-${_MKL_VER}/mkl") # default until ParallelStudioXE2015 (no root permissions)
	endforeach()
    if(APPLE)
	    list(APPEND _mkl_search_paths "/opt/intel/compilers_and_libraries/mac/mkl") # default for ParallelStudioXE2016 and later (root permissions)
    else()
        list(APPEND _mkl_search_paths "/opt/intel/compilers_and_libraries/linux/mkl") # default for ParallelStudioXE2016 and later (root permissions)
    endif()
	list(APPEND _mkl_search_paths "/opt/intel/oneapi/mkl") # default for oneAPI (2020) and later (root permissions)
    if(APPLE)
	    list(APPEND _mkl_search_paths "$ENV{HOME}/intel/compilers_and_libraries/mac/mkl") # default for ParallelStudioXE2016 and later (no root permissions)
    else()
        list(APPEND _mkl_search_paths "$ENV{HOME}/intel/compilers_and_libraries/linux/mkl") # default for ParallelStudioXE2016 and later (no root permissions)
    endif()
    list(APPEND _mkl_search_paths "$ENV{HOME}/intel/oneapi/mkl") # default for oneAPI (2020) and later (no root permissions)
endif()

if(CMAKE_MKL_DEBUG)
    message(STATUS "Searching for MKL in:")
    foreach (dir ${_mkl_search_paths})
        message(${dir})
    endforeach ()         
endif(CMAKE_MKL_DEBUG)

# Determine MKL's library folder
#
set(_mkl_libpath_suffix_orig "lib")

set(_mkl_libpath_suffix  ${_mkl_libpath_suffix_orig})
list(APPEND _mkl_libpath_suffix "latest/${_mkl_libpath_suffix_orig}")

if(WIN32)
    set(_mkl_platform_suffix "_win")
    set(_mkl_libname_prefix "")
    set(_mkl_shared_lib "_dll.lib")
    set(_mkl_static_lib ".lib")
elseif(APPLE)
    set(_mkl_platform_suffix "_mac")
    set(_mkl_libname_prefix "lib")
    set(_mkl_shared_lib ".dylib")
    set(_mkl_static_lib ".a")
else() # LINUX
    set(_mkl_platform_suffix "_lin")
    set(_mkl_libname_prefix "lib")
    set(_mkl_shared_lib ".so")
    set(_mkl_static_lib ".a")
endif()
if(CMAKE_SIZEOF_VOID_P EQUAL 4) # 32 bit
    set(_mkl_arch_suffix "ia32")
else() # 64 bit
    set(_mkl_arch_suffix "intel64")
endif()
list(APPEND _mkl_libpath_suffix "${_mkl_libpath_suffix_orig}/${_mkl_arch_suffix}")
list(APPEND _mkl_libpath_suffix "latest/${_mkl_libpath_suffix_orig}/${_mkl_arch_suffix}")
list(APPEND _mkl_libpath_suffix "${_mkl_libpath_suffix_orig}/${_mkl_arch_suffix}${_mkl_plaform_suffix}")
list(APPEND _mkl_libpath_suffix "latest/${_mkl_libpath_suffix_orig}/${_mkl_arch_suffix}${_mkl_plaform_suffix}")

foreach (dir IN LISTS ${EXTRA_SEARCH_DIRS})
  list(APPEND _mkl_search_paths "${dir}" "${dir}/mkl" "${dir}/mkl/compiler")
endforeach ()              

if(CMAKE_MKL_DEBUG)
message(STATUS "MKL_LIBPATH_SUFFIX: ${_mkl_libpath_suffix}")
endif(CMAKE_MKL_DEBUG)
# Functions: finds both static and shared MKL libraries
#
function(__mkl_find_library _varname _libname)
    find_library(${_varname}_DYN
          NAMES ${_mkl_libname_prefix}${_libname}${_mkl_shared_lib}
          HINTS ${_mkl_search_paths}
          PATH_SUFFIXES ${_mkl_libpath_suffix})
    mark_as_advanced(${_varname}_DYN)
    find_library(${_varname}_ST
          NAMES ${_mkl_libname_prefix}${_libname}${_mkl_static_lib}
          HINTS ${_mkl_search_paths}
          PATH_SUFFIXES ${_mkl_libpath_suffix})
    mark_as_advanced(${_varname}_ST)
endfunction()

# Find MKL headers
#
find_path(MKL_INCLUDE_DIR mkl.h
    PATHS ${_mkl_search_paths}
    PATH_SUFFIXES
    "include" #older versions
    "latest/include" #newer (2021) oneAPI version
    "include/mkl" #when installed by apt-get on debian
    )
mark_as_advanced(MKL_INCLUDE_DIR)

if(CMAKE_MKL_DEBUG)  
  message(STATUS "MKL INCLUDE DIR ${MKL_INCLUDE_DIR}")
endif()
# Group flags for static libraries on Linux (GNU, PGI, ICC -> same linker)
#
if(UNIX AND NOT APPLE)
    set(_mkl_linker_pre_flags "-Wl,--start-group")
    set(_mkl_linker_post_flags "-Wl,--end-group -ldl -lm")
endif()

# Core MKL
#
__mkl_find_library(MKL_CORE_LIB mkl_core)
if(CMAKE_MKL_DEBUG)  
  message(STATUS "MKL CORE LIB ${MKL_CORE_LIB_DYN}")
endif()
# Interface
#
__mkl_find_library(MKL_INTERFACE_INTEL_32BIT_LIB mkl_intel_lp64)
__mkl_find_library(MKL_INTERFACE_INTEL_64BIT_LIB mkl_intel_ilp64)
if(NOT APPLE AND CMAKE_Fortran_COMPILER_LOADED
             AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    __mkl_find_library(MKL_INTERFACE_GF_32BIT_LIB mkl_gf_lp64)
    __mkl_find_library(MKL_INTERFACE_GF_64BIT_LIB mkl_gf_ilp64)
endif()
# Threading
#
__mkl_find_library(MKL_SEQ_LIB mkl_sequential)
if(NOT APPLE AND (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR
                  CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
                  CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"))
    __mkl_find_library(MKL_OMP_LIB mkl_gnu_thread)
else()
    __mkl_find_library(MKL_OMP_LIB mkl_intel_thread)
endif()
__mkl_find_library(MKL_TBB_LIB mkl_tbb_thread)
if(CMAKE_MKL_DEBUG)
  message(STATUS "MKL THREAD LIB ${MKL_OMP_LIB_DYN}")
endif()
# BLACS
#
if(APPLE)
    __mkl_find_library(MKL_BLACS_MPICH_32BIT_LIB mkl_blacs_mpich_lp64)
    __mkl_find_library(MKL_BLACS_MPICH_64BIT_LIB mkl_blacs_mpich_ilp64)
else()
    __mkl_find_library(MKL_BLACS_MPICH_32BIT_LIB mkl_blacs_intelmpi_lp64)
    __mkl_find_library(MKL_BLACS_MPICH_64BIT_LIB mkl_blacs_intelmpi_ilp64)
endif()
__mkl_find_library(MKL_BLACS_OMPI_32BIT_LIB mkl_blacs_openmpi_lp64)
__mkl_find_library(MKL_BLACS_OMPI_64BIT_LIB mkl_blacs_openmpi_ilp64)

# ScaLAPACK
#
__mkl_find_library(MKL_SCALAPACK_32BIT_LIB mkl_scalapack_lp64)
__mkl_find_library(MKL_SCALAPACK_64BIT_LIB mkl_scalapack_ilp64)

# Check if core libs were found
#
find_package_handle_standard_args(MKL REQUIRED_VARS MKL_INCLUDE_DIR
                                                    Threads_FOUND)

# Sequential has no threading dependency. There is currently no TBB module
# shipped with CMake. The dependency is not accounted for.
#
set(_mkl_dep_found_SEQ TRUE)
if(TARGET TBB::tbb)
  set(_mkl_dep_TBB TBB::tbb)
  set(_mkl_dep_found_TBB TRUE)
endif()
if (TARGET OpenMP::OpenMP_CXX)
  set(_mkl_dep_OMP OpenMP::OpenMP_CXX)
  set(_mkl_dep_found_OMP TRUE)
endif()

# Define all blas, blacs and scalapack
#
foreach(_libtype "ST" "DYN")
    set(_mkl_core_lib ${MKL_CORE_LIB_${_libtype}})
    foreach(_bits "32BIT" "64BIT")
        set(_mkl_scalapack_lib ${MKL_SCALAPACK_${_bits}_LIB_${_libtype}})
        foreach(_iface "INTEL" "GF")
            set(_mkl_interface_lib ${MKL_INTERFACE_${_iface}_${_bits}_LIB_${_libtype}})
            foreach(_threading "SEQ" "OMP" "TBB")
                set(_mkl_threading_lib ${MKL_${_threading}_LIB_${_libtype}})

                string(TOLOWER "${_iface}_${_bits}_${_threading}_${_libtype}" _tgt_config)
                set(_mkl_tgt mkl::mkl_${_tgt_config})

                if(MKL_FOUND
                   AND _mkl_interface_lib
                   AND _mkl_threading_lib
                   AND _mkl_core_lib
                   AND _mkl_dep_found_${_threading}
                   AND NOT TARGET ${_mkl_tgt})
                    set(_mkl_libs "${_mkl_linker_pre_flags}"
                                  "${_mkl_interface_lib}"
                                  "${_mkl_threading_lib}"
                                  "${_mkl_core_lib}"
                                  "${_mkl_dep_${_threading}}"
                                  "${_mkl_linker_post_flags}"
                                  "Threads::Threads")
                    add_library(${_mkl_tgt} INTERFACE IMPORTED)
                    set_target_properties(${_mkl_tgt} PROPERTIES
                      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}"
                      INTERFACE_LINK_LIBRARIES "${_mkl_libs}")
                    if(CMAKE_MKL_DEBUG)
                        message(STATUS "MKL_TGT: ${_mkl_tgt}")
                        message(STATUS "_mkl_libs: ${_mkl_libs}")
                        message(STATUS "_mkl_dep_${_threading}: ${_mkl_dep_${_threading}}")

                    endif(CMAKE_MKL_DEBUG)
                endif()

                foreach(_mpi_impl "MPICH" "OMPI")
                    set(_mkl_blacs_lib ${MKL_BLACS_${_mpi_impl}_${_bits}_LIB_${_libtype}})

                    string(TOLOWER "${_mpi_impl}_${_iface}_${_bits}_${_threading}_${_libtype}" _tgt_config)
                    set(_blacs_tgt mkl::blacs_${_tgt_config})
                    set(_scalapack_tgt mkl::scalapack_${_tgt_config})

                    if(_mkl_blacs_lib
                       AND TARGET ${_mkl_tgt}
                       AND TARGET MPI::MPI_CXX
                       AND NOT TARGET ${_blacs_tgt})
                        set(_blacs_libs "${_mkl_linker_pre_flags}"
                                        "${_mkl_interface_lib}"
                                        "${_mkl_threading_lib}"
                                        "${_mkl_core_lib}"
                                        "${_mkl_blacs_lib}"
                                        "${_mkl_linker_post_flags}"
                                        "MPI::MPI_CXX"
                                        "${_mkl_dep_${_threading}}"
                                        "Threads::Threads")
                        add_library(${_blacs_tgt} INTERFACE IMPORTED)
                        set_target_properties(${_blacs_tgt} PROPERTIES
                            INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}"
                            INTERFACE_LINK_LIBRARIES "${_blacs_libs}")
                    endif()

                    if(_mkl_scalapack_lib
                       AND TARGET ${_blacs_tgt}
                       AND NOT TARGET ${_scalapack_tgt})
                        set(_scalapack_libs "${_mkl_scalapack_lib}"
                                            "${_blacs_tgt}")
                        add_library(${_scalapack_tgt} INTERFACE IMPORTED)
                        set_target_properties(${_scalapack_tgt} PROPERTIES
                            INTERFACE_LINK_LIBRARIES "${_scalapack_libs}")
                    endif()
                endforeach()
            endforeach()
        endforeach()
    endforeach()
endforeach()
mark_as_advanced(MKL_ROOT)
