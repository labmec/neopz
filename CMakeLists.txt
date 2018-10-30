# @file neopz/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.11.0)

#disable excessive warning messages by Microsoft VC compiler
set (CMAKE_USER_MAKE_RULES_OVERRIDE "${CMAKE_CURRENT_LIST_DIR}/CompilerOptionsMSVC.cmake")

set(CMAKE_CONFIGURATION_TYPES Debug Release CACHE TYPE  INTERNAL FORCE )

project (PZ)

#set(CMAKE_POSITION_INDEPENDENT_CODE ON)

IF (WIN32)
	# Prevents timespec redefinition problem with Visual Studio 2015
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
ELSE (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
	if ( NOT APPLE AND NOT CMAKE_BUILD_TYPE )
		set (CMAKE_BUILD_TYPE Release CACHE STRING "Options: Release Debug None MinSizeRel RelWithDebInfo" FORCE) 
	endif()
ENDIF()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## inserting main definitions
#specifying the path to neopz source code
set(PZSOURCEDIR ${PROJECT_SOURCE_DIR})
#specify where the refinement patterns can be found (default behaviour when running the code)
set(REFPATTERNDIR ${PROJECT_SOURCE_DIR}/Refine/RefPatterns)


#Defining variables
#To cmake
string(FIND ${CMAKE_CXX_FLAGS_DEBUG} "-DPZDEBUG" val)
if(${val} LESS 0)
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DPZDEBUG" CACHE STRING "c++ debug flags go here" FORCE)
endif()

string(FIND ${CMAKE_CXX_FLAGS_RELEASE} "-DNODEBUG" val)
if(${val} LESS 0)
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNODEBUG -DZERO_INTERNAL_RESIDU -DMAKEINTERNAL" CACHE STRING "c++ release flags go here" FORCE)
endif()

# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/local/opt/qt5/lib/cmake/ /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/pthread ./../externallibs/metis ./../externallibs/metis/include ./../externallibs/metis/lib ./..metis ./../metis/include ./../metis/lib ./externallibs/lib ./externallibs/include ./../externallibs/lib ./../externallibs/include ./../externallibs/include/lapack /usr/lib/atlas-base/ /usr/include/atlas /usr/local/atlas/include /usr/local/atlas/lib )

# INSTALL part
#To install neopz/include
set(NEOPZ_INCLUDE pzlib/include)
#To install neopz library
set(NEOPZ_LIB pzlib/lib)

#ENABLING EXTERNAL LIBRARIES

#enabling pthread library
find_path(PTHREAD_INCLUDE NAMES PThread.h pthread.h PATHS ${SEARCH_DIRS} ${PROJECT_SOURCE_DIR}/../externallibs/pthread ${PROJECT_SOURCE_DIR}/../externallibs/pthread/include ${PROJECT_SOURCE_DIR}/../pthread )
find_library(PTHREAD_LIB NAMES libpthread.dylib pthread.lib pthreadVC2.lib libPThread.so libPThread.a libpthread.so libpthread.a PATHS ${SEARCH_DIRS} ${PROJECT_SOURCE_DIR}/../externallibs/pthread/lib  ${PROJECT_SOURCE_DIR}/../pthread/lib)
if(PTHREAD_INCLUDE-NOTFOUND)
	set (PTHREAD_INCLUDE "" CACHE PATH "Directory where pthread.h can be found")
endif()
if(PTHREAD_LIB-NOTFOUND)
	set (PTHREAD_LIB "" CACHE PATH "Directory where the pthread library can be found")
endif()
include_directories(${PTHREAD_INCLUDE})

IF(APPLE)
   # Mac OS X specific code
   add_definitions(-DMACOSX)
   set(MACOSX MACOSX)
   if (NOT CMAKE_GENERATOR MATCHES "Ninja")
       FIND_LIBRARY(ACCELERATE_LIB Accelerate)
	   link_libraries(${ACCELERATE_LIB})
   else()
       set(CMAKE_FIND_FRAMEWORK NEVER)
	   find_library(ACCELERATE_LIB Accelerate PATHS /System/Library/Frameworks/Accelerate.framework/Accelerate)
	   set(CMAKE_FIND_FRAMEWORK FIRST)
   endif()
ENDIF(APPLE)

#enabling Boost library
option(USING_BOOST "Whether the Boost package will be linked with the library" OFF)

#To enable unit test using boost UnitTest framework
option(BUILD_UNITTESTING "Whether the Boost UnitTest will be linked in" OFF)
if(BUILD_UNITTESTING)
	set(USING_BOOST ON CACHE BOOL "Whether the Boost package will be linked with the library" FORCE)
endif()

if(USING_BOOST)
	IF (WIN32)
		set(Boost_USE_STATIC_LIBS   ON)
		set(Boost_USE_MULTITHREADED ON)
		set(BOOST_ROOT "../externallibs/lib/boost_1_65_0")	
	ENDIF()
	find_package( Boost COMPONENTS unit_test_framework date_time atomic REQUIRED )
endif(USING_BOOST)

if(BUILD_UNITTESTING)
	IF (NOT WIN32)
		link_libraries(${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
	ENDIF()
endif()

option(USING_FAD "Whether the FAD package will be linked with the library" OFF)
if(USING_FAD)
	#Adding a variable to hold this definition
	set(_AUTODIFF _AUTODIFF)
endif(USING_FAD)

#enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
if(USING_LOG4CXX)
	#Adding a variable to hold this definition
	set(LOG4CXX LOG4CXX)
	find_path(LOG4CXX_INCLUDE log4cxx PATHS ${SEARCH_DIRS} ${PROJECT_SOURCE_DIR}/../log4cxx ${PROJECT_SOURCE_DIR}/../log4cxx_src/include ${PROJECT_SOURCE_DIR}/../externallibs/log4cxx_src ${PROJECT_SOURCE_DIR}/../externallibs/log4cxx_src/include)
	find_library(LOG4CXX_LIB_RELEASE
		         NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a
				 PATHS ${SEARCH_DIRS} ./../log4cxx_src ./../log4cxx_src/lib ./../externallibs/log4cxx_src ./../externallibs/log4cxx_src/lib ./../externallibs/lib ./../externallibs/lib/Release 
				)
	find_library(LOG4CXX_LIB_DEBUG
				 NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a       # Change name here if debug lib has different name.
				 PATHS ${SEARCH_DIRS} ./../log4cxx_src ./../log4cxx_src/lib ./../externallibs/log4cxx_src ./../externallibs/log4cxx_src/lib ./../externallibs/lib ./../externallibs/lib/Debug
				)
	set(LOG4CXX_LIB
	    debug ${LOG4CXX_LIB_DEBUG}
	    optimized ${LOG4CXX_LIB_RELEASE}
	   )	
	if(LOG4CXX_INCLUDE-NOTFOUND)
		set (LOG4CXX_INCLUDE "" CACHE PATH "Directory where log4cxx header files can be found")
	else()
		include_directories(${LOG4CXX_INCLUDE})
	endif()
	
	if(LOG4CXX_LIB-NOTFOUND)
		set (LOG4CXX_LIB "" CACHE PATH "Directory where the log4cxx library can be found")
	endif()
endif(USING_LOG4CXX)

#enabling compiler interprocedural optimizations (IPO)
option(USING_IPO "Whether link and assembly support for compiler IPO will be activated or not (may cause linking to be very slow)." OFF)

#enabling METIS library
option(USING_METIS "Whether the METIS library will be linked in" OFF)
if(USING_METIS)
	####add_definitions(-DUSING_METIS)
	find_path(METIS_INCLUDE metis.h PATHS ${SEARCH_DIRS})
	find_library(METIS_LIB NAMES libmetis.dylib metis.lib libmetis.so libmetis.a PATHS ${SEARCH_DIRS})
	if(METIS_INCLUDE-NOTFOUND)
		set (METIS_INCLUDE "" CACHE PATH "Directory where metis header files can be found")
	else()
		include_directories(${METIS_INCLUDE})
	endif()
	if(METIS_LIB-NOTFOUND)
		set (METIS_LIB "" CACHE PATH "Directory where the metis library can be found")
	endif()
endif(USING_METIS)

#enabling OPENSSL library (used for MD5 computation)
option(USING_OPENSSL "Whether the OPENSSL library will be linked in" OFF)
if(USING_OPENSSL)
	link_libraries(crypto)
endif(USING_OPENSSL)

#enabling TBB library
option(USING_TBB "Whether the TBB library will be linked in" OFF)
if(USING_TBB)
	if(WIN32)
		set(TBB_ROOT "C:/Program Files (x86)/IntelSWTools/TBB")
		find_path(TBB_INCLUDE tbb.h HINTS ${TBB_ROOT}/tbb2018_20180618oss/include/tbb NO_DEFAULT_PATH)		
		find_library(TBB_LIB_RELEASE tbb.lib HINTS ${TBB_ROOT}/tbb2018_20180618oss/lib/intel64/vc14 NO_DEFAULT_PATH)
		find_library(TBB_LIB_DEBUG tbb_debug.lib HINTS ${TBB_ROOT}/tbb2018_20180618oss/lib/intel64/vc14 NO_DEFAULT_PATH)
		
		set(TBB_LIB
      		debug ${TBB_LIB_DEBUG}
	        optimized ${TBB_LIB_RELEASE}
	    )	
		
		if(TBB_INCLUDE-NOTFOUND)
			set (TBB_INCLUDE "" CACHE PATH "Directory where tbb header files can be found")
		else()
			include_directories(${TBB_INCLUDE}/..)
		endif()
		if(TBB_LIB-NOTFOUND)
			set (TBB_LIB "" CACHE PATH "Directory where the tbb library can be found")
		else()
			link_libraries(${TBB_LIB})
		endif()
	else()
		find_file(TBB_HEADER_FILE tbb.h HINTS /usr/local/TBB/*/include/tbb/ NO_DEFAULT_PATH)

		if(NOT TBB_HEADER_FILE MATCHES "tbb.h")
		  set (TBB_HEADER_FILE "" CACHE PATH "Path to 'tbb.h' file")
		  message( FATAL_ERROR "TBB_HEADER_FILE must contain a valid path to 'tbb.h' file!")
		else()
		  get_filename_component(TBB_ROOT_DIR ${TBB_HEADER_FILE} DIRECTORY)
		  get_filename_component(TBB_INCLUDE ${TBB_ROOT_DIR} DIRECTORY)
		  get_filename_component(TBB_ROOT_DIR ${TBB_INCLUDE} DIRECTORY)
		  if(APPLE)
			set(TBB_LIB ${TBB_ROOT_DIR}/lib/libtbb.dylib)
		  else() # Linux
			set(TBB_LIB ${TBB_ROOT_DIR}/lib/intel64/gcc4.7/libtbb.so)
		  endif()

		  include_directories(${TBB_INCLUDE})
		  link_libraries(${TBB_LIB})
		endif()
	endif(WIN32)
endif(USING_TBB)

#enabling OpenMP runtime and library
option(USING_OPENMP "Whether the OpenMP library and runtime will be linked in" OFF)
if(USING_OPENMP)
	find_path(OPENMP_INCLUDE omp.h PATHS ${SEARCH_DIRS} /opt/local/lib/gcc5/gcc/x86_64-apple-darwin15/5.2.0/include)
	find_library(OPENMP_LIB NAMES libgomp.dylib PATHS ${SEARCH_DIRS}  /opt/local/lib/gcc5)
	if(OPENMP_INCLUDE-NOTFOUND)
		set (OPENMP_INCLUDE "" CACHE PATH "Directory where omp header files can be found")
	else()
		include_directories(${OPENMP_INCLUDE})
	endif()
	if(OPENMP_LIB-NOTFOUND)
		set (TBB_OPENMP "" CACHE PATH "Directory where the OpenMP library can be found")
	else()
		link_libraries(${OPENMP_LIB})
	endif()
endif(USING_OPENMP)

#enabling LIKWID library
option(USING_LIKWID "Whether the likwid library will be linked in" OFF)
if(USING_LIKWID)
	find_path(LIKWID_INCLUDE likwid.h PATHS ${SEARCH_DIRS})
	find_library(LIKWID_LIB NAMES liblikwid.a PATHS ${SEARCH_DIRS})
	if(LIKWID_INCLUDE-NOTFOUND)
		set (LIKWID_INCLUDE "" CACHE PATH "Directory where likwid header files can be found")
	else()
		include_directories(${LIKWID_INCLUDE})
	endif()
	if(LIKWID_LIB-NOTFOUND)
		set (LIKWID_LIB "" CACHE PATH "Directory where the likwid library can be found")
	else()
		link_libraries(${LIKWID_LIB})
	endif()
endif(USING_LIKWID)

#enabling numa library
option(USING_LIBNUMA "Whether the Lib NUMA library will be linked in" OFF)
if(USING_LIBNUMA)
	find_path(LIBNUMA_INCLUDE numa.h PATHS ${SEARCH_DIRS})
	find_library(LIBNUMA_LIB NAMES libnuma.so PATHS ${SEARCH_DIRS})
	if(LIBNUMA_INCLUDE-NOTFOUND)
		set (LIBNUMA_INCLUDE "" CACHE PATH "Directory where libnuma header files can be found (consider installing libnuma-dev)")
	else()
		include_directories(${LIBNUMA_INCLUDE})
	endif()
	if(LIBNUMA_LIB-NOTFOUND)
		set (LIBNUMA_LIB "" CACHE PATH "Directory where the libnuma library can be found")
	else()
		link_libraries(${LIBNUMA_LIB})
	endif()
endif(USING_LIBNUMA)

# ----------------------- enabling MatLab library -----------------------
option(USING_MATLAB_ENGINE "Whether MatLab Engine will be linked in" OFF)
if(USING_MATLAB_ENGINE)
	find_package(MATLAB REQUIRED)
	if (NOT MATLAB_FOUND)
		message(FATAL_ERROR "This project requires the MATLAB library. Please inform the path to your matlab application")
	endif()

	if(NOT EXISTS "/usr/bin/matlab")
		message(FATAL_ERROR "Error! Please create matlab symlink on /usr/bin by typing on terminal:\nsudo ln -sf PATH_TO_YOUR_MATLAB/bin/matlab /usr/bin/matlab")
	endif()

	include_directories(${Matlab_ROOT_DIR}/extern/include)

endif(USING_MATLAB_ENGINE)

#enabling LAPACK library
option(USING_LAPACK "Whether the LAPACK library will be linked in" OFF)
option(USING_BLAS "Whether the Basic Linear Algebra Subprograms (BLAS) library will be linked in" OFF)
if(USING_LAPACK)
	set(USING_BLAS ON CACHE BOOL "Whether the Basic Linear Algebra Subprograms (BLAS) library will be linked in" FORCE)
	if(APPLE)
	else()
		find_path(LAPACK_INCLUDE lapacke.h PATHS ${SEARCH_DIRS})
		find_library(LAPACK_LIB NAMEs liblapack.lib liblapack.so PATHS ${SEARCH_DIRS})
		if(LAPACK_INCLUDE-NOTFOUND)
			set (LAPACK_INCLUDE "" CACHE PATH "Directory where LAPACK header files can be found")
		else()
			#include_directories(${LAPACK_INCLUDE})
		endif()
		if(LAPACK_LIB-NOTFOUND)
			set (LAPACK_LIB "" CACHE PATH "Directory where the LAPACK library can be found")
		endif()
	endif(APPLE)
endif(USING_LAPACK)

#enabling BLAS library
#option(USING_BLAS "Whether the Basic Linear Algebra Subprograms (BLAS) library will be linked in" OFF)
if(USING_BLAS)
	if(APPLE)
	else()
		find_path(CBLAS_INCLUDE cblas.h PATHS ${SEARCH_DIRS})
		find_library(CBLAS_LIB NAMES libcblas cblas libtcblas PATHS ${SEARCH_DIRS})
		if(CBLAS_INCLUDE-NOTFOUND)
			set (CBLAS_INCLUDE "" CACHE PATH "Directory where CBLAS header files can be found")
		else()
			include_directories(${CBLAS_INCLUDE})
		endif()
		if(CBLAS_LIB-NOTFOUND)
			set (CBLAS_LIB "" CACHE PATH "Directory where the CBLAS library can be found")
		endif()
	endif(APPLE)
endif(USING_BLAS)


#enabling the new skyline matrix format
option(USING_NEW_SKYLMAT "Enable the new skyline matrix format" OFF)
if(USING_NEW_SKYLMAT)
		set(USING_NEW_SKYLMAT DUSING_NEW_SKYLMAT)
endif(USING_NEW_SKYLMAT)

#enabling PAPI library
option(USING_PAPI "Whether the PAPI library will be linked in" OFF)
if(USING_PAPI)
        find_path(PAPI_INCLUDE NAMES papi.h PATHS ${SEARCH_DIRS})
        find_library(PAPI_LIB NAMES libpapi.a libpapi.so PATHS ${SEARCH_DIRS})
        if(PAPI_INCLUDE-NOTFOUND)
                set (PAPI_INCLUDE "" CACHE PATH "Directory where papi.h can be found")
        endif()
        if(PAPI_LIB-NOTFOUND)
                set (PAPI_LIB "" CACHE PATH "Directory where the papi library can be found")
        endif()
endif(USING_PAPI)

#enabling HWLOC library
option(USING_HWLOC "Whether the HWLOC library will be linked in" OFF)
if(USING_HWLOC)
        find_path(HWLOC_INCLUDE NAMES hwloc.h PATHS ${SEARCH_DIRS})
        find_library(HWLOC_LIB NAMES hwloc.a libhwloc.so libhwloc.dylib PATHS ${SEARCH_DIRS})
        if(HWLOC_INCLUDE-NOTFOUND)
                set (HWLOC_INCLUDE "" CACHE PATH "Directory where hwloc.h can be found")
        endif()
        if(HWLOC_LIB-NOTFOUND)
                set (HWLOC_LIB "" CACHE PATH "Directory where the hwloc library can be found")
        endif()
endif(USING_HWLOC)



#enabling MKL library
option(USING_MKL "Whether the MKL library will be linked in" OFF)
if(USING_MKL)
	if (WIN32)
		set(MKL_ROOT "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows")
		find_path(MKL_INCLUDE NAMES mkl.h PATHS ${MKL_ROOT}/mkl/include)
		find_library(MKL_LIB_INTEL NAMES mkl_rt.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
		find_library(MKL_LIB_CORE NAMES mkl_core.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
		find_library(MKL_LIB_THREAD NAMES mkl_intel_thread.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
                find_library(COMPOSER_OMP NAMES libiomp5md.lib PATHS ${MKL_ROOT}/compiler/lib/intel64_win)
	else()
		find_path(MKL_INCLUDE NAMES mkl.h PATHS ${SEARCH_DIRS} /opt/intel/mkl/include /softwares/intel/mkl/include)
		find_library(MKL_LIB_INTEL NAMES libmkl_intel_lp64.so libmkl_intel_lp64.dylib PATHS
				${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64)
		find_library(MKL_LIB_CORE NAMES libmkl_core.so libmkl_core.dylib PATHS ${SEARCH_DIRS}
				/opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
		find_library(MKL_LIB_THREAD NAMES libmkl_intel_thread.so libmkl_intel_thread.dylib PATHS
				${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
		find_library(COMPOSER_OMP NAMES libiomp5.so libiomp5.dylib PATHS ${SEARCH_DIRS}
				/opt/intel/composer_xe/compiler/lib
				/opt/intel/composerxe/lib/intel64
				/opt/intel/lib /softwares/intel/lib/intel64
                                /opt/intel/compilers_and_libraries/linux/lib/intel64_lin )
	endif()

	if(MKL_INCLUDE-NOTFOUND)
                set (MKL_INCLUDE "" CACHE PATH "Directory where mkl.h can be found")
        else()
                #include_directories(${MKL_INCLUDE})
        endif()
        if(MKL_LIB-NOTFOUND)
                set (MKL_LIB_INTEL "" CACHE PATH "Directory where the mkl library can be found")
        else()
                link_libraries(${MKL_LIB_INTEL})
                link_libraries(${MKL_LIB_CORE})
                link_libraries(${MKL_LIB_THREAD})
                link_libraries(${COMPOSER_OMP})
        endif()
endif(USING_MKL)

#To enable Python bindings
option(BUILD_PYTHON_BINDING "To build the Python bindings" OFF)
if (BUILD_PYTHON_BINDING)
  find_package(SWIG "2.0")
  find_package(PythonLibs)
  if (NOT SWIG_FOUND OR NOT PYTHONLIBS_FOUND)
    message(WARNING "SWIG 2 and Python headers required for building Python API bindings. You can set proper SWIG executable with argument -DSWIG_EXECUTABLE=/path/to/swig2.0")
    set(BUILD_PYTHON_BINDING OFF)
  endif()
endif()

#To building several examples for training and experimentation PZ
option(BUILD_TUTORIAL "To build examples for training or teaching PZ" OFF)
if(BUILD_TUTORIAL)
	####add_definitions(-DBUILD_TUTORIAL)
endif()

#To building several examples for training and experimentation PZ
option(BUILD_BROKEN "To build projects which are broken" OFF)

if(BUILD_TUTORIAL)
	####add_definitions(-DBUILD_TUTORIAL)
endif()


#To building all the specific projects developed with PZ environment
option(BUILD_PROJECTS "To build specific projects." OFF)

#To building all the performance projects developed with PZ environment
option(BUILD_PERF_TESTS "To build the performance test." OFF)

#To building all the plasticity materials
option(BUILD_PLASTICITY_MATERIALS "To build plasticity materials." OFF)

#To select type of REAL
set(REALTypesAre "Types are: double, float, long double, pzfpcounter")
set(REAL_TYPE "double" CACHE STRING "${REALTypesAre}")

if (REAL_TYPE STREQUAL "double")
  set(REAL_TYPE_DEF "REALdouble")
elseif (REAL_TYPE STREQUAL "float")
  set(REAL_TYPE_DEF "REALfloat")
elseif (REAL_TYPE STREQUAL "long double")
  set(REAL_TYPE_DEF "REALlongdouble")
elseif (REAL_TYPE STREQUAL "pzfpcounter")
  set(REAL_TYPE_DEF "REALpzfpcounter")
else()
  message (FATAL_ERROR "ERROR: Please specify a valid type for REAL. ${REALTypesAre}")
endif()
message (STATUS "Using REAL type = ${REAL_TYPE}")

#To select type of STATE
set(STATETypesAre "Types are: double, float, long double, complex<float>, complex<double>, complex<long double>")
set(STATE_TYPE "double" CACHE STRING "${STATETypesAre}")

if (STATE_TYPE STREQUAL "double")
  set(STATE_TYPE_DEF "STATEdouble")
elseif (STATE_TYPE STREQUAL "float")
  set(STATE_TYPE_DEF "STATEfloat")
elseif (STATE_TYPE STREQUAL "long double")
  set(STATE_TYPE_DEF "STATElongdouble")
elseif (STATE_TYPE STREQUAL "complex<float>")
  set(STATE_TYPE_DEF "STATEcomplexf")
  set(STATE_COMPLEX "STATE_COMPLEX")
  set (BUILD_COMPLEX_PROJECTS ON)
elseif (STATE_TYPE STREQUAL "complex<double>")
  set(STATE_TYPE_DEF "STATEcomplexd")
  set(STATE_COMPLEX "STATE_COMPLEX")
  set (BUILD_COMPLEX_PROJECTS ON)
elseif (STATE_TYPE STREQUAL "complex<long double>")
  set(STATE_TYPE_DEF "STATEcomplexld")
  set(STATE_COMPLEX "STATE_COMPLEX")
  set (BUILD_COMPLEX_PROJECTS ON)
else()
  message (FATAL_ERROR "ERROR: Please specify a valid type for STATE. ${STATETypesAre}")
endif()
message (STATUS "Using STATE type = ${STATE_TYPE}")

# Some definitions for windows
if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on pz_config.h")
	add_definitions(-DVC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	#disabling STL warning for std::_Vector_alloc when LOG4CXX is enabled
	#if (USING_LOG4CXX)
	    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /wd4251 /wd4244")
	    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /wd4251 /wd4244")
	#endif()
	#define use of pthread static lib.
	add_definitions(-DPTW32_STATIC_LIB)
    #enabling /bigobj
    add_definitions("/bigobj")
	if (NOT CMAKE_GENERATOR MATCHES "Ninja")
		add_definitions("/MP")
    endif()
    FOREACH(FLAG_TYPE EXE MODULE SHARED)
		# preventing .pdb unnecessary file generation overhead: no /DEBUG flag in linker flags
        SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "/INCREMENTAL:NO" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
        MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
    ENDFOREACH ()

endif()

if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
  find_package(Git)
  if(GIT_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "PZ_BRANCH"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "PZ_REVISION"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} log -1 --format=%cd --date=local
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
      OUTPUT_VARIABLE "PZ_REVISION_DATE"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  else()
    set(PZ_BRANCH "unknown")
	set(PZ_REVISION "unknown")
	set(PZ_REVISION_DATE "unknown")
  endif(GIT_FOUND)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)

# Begin of pz_config.h part
# This file contains all definitions used to compile the library pz
# this will be saved at include diretory in the install target
set (CONFIG_H "${PROJECT_BINARY_DIR}/Common/pz_config.h")
configure_file ("${PROJECT_SOURCE_DIR}/Common/pz_config.h.in" ${CONFIG_H})
# add it to include search path
include_directories("${PROJECT_BINARY_DIR}/Common")
# add the .h file in the install target
install (FILES ${CONFIG_H} DESTINATION ${NEOPZ_INCLUDE} COMPONENT libraries)
# End of pz_config.h part


#creating target pz library and installing it
if (BUILD_PYTHON_BINDING)
	IF (NOT WIN32)
		#Build a shared library to be used with python binding
		add_library(pz SHARED)
	ELSE()
		add_library(pz STATIC)
	ENDIF()
else()
  #Build a static library (default)
  add_library(pz)
endif()


target_link_libraries(pz ${PTHREAD_LIB})
target_include_directories(pz PUBLIC ${PTHREAD_INCLUDE})

IF (USING_IPO)
  set_property(DIRECTORY PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
ENDIF()

IF (NOT WIN32)
	target_link_libraries(pz pthread)
ENDIF()

if (USING_BOOST)
   target_link_libraries(pz ${Boost_LIBRARIES})
   target_include_directories(pz PUBLIC ${Boost_INCLUDE_DIR})
endif()

if (USING_LOG4CXX)
  target_link_libraries(pz ${LOG4CXX_LIB})
  target_include_directories(pz PUBLIC ${LOG4CXX_INCLUDE})
  if(WIN32)
	target_link_libraries(pz odbc32.lib ws2_32.lib mswsock.lib)
  endif()
endif()

if (USING_OPENMP)
  target_link_libraries(pz ${OPENMP_LIB})
  target_include_directories(pz PUBLIC ${OPENMP_INCLUDE})
endif()

if (USING_METIS)
  target_link_libraries(pz ${METIS_LIB})
  target_include_directories(pz PUBLIC ${METIS_INCLUDE})
endif()

if (USING_LAPACK)
  target_link_libraries(pz ${LAPACK_LIB} )
  target_include_directories(pz PUBLIC ${LAPACK_INCLUDE})
endif()

if (USING_BLAS)
  target_link_libraries(pz ${CBLAS_LIB} ) # lapack BLAS
endif()

if (USING_PAPI)
  target_link_libraries(pz ${PAPI_LIB})
endif()

if (USING_MKL)
  target_link_libraries(pz ${MKL_LIB_INTEL} ${MKL_LIB_CORE} ${MKL_LIB_THREAD} ${COMPOSER_OMP})
  target_include_directories(pz PUBLIC ${MKL_INCLUDE})
endif()

if (USING_TBB)
  target_link_libraries(pz ${TBB_LIB} )
  target_include_directories(pz PUBLIC ${TBB_INCLUDE})
endif()


if (USING_HWLOC)
  target_link_libraries(pz ${HWLOC_LIB})
endif()

#increment all the directories where the sources files is
if(NOT BUILD_COMPLEX_PROJECTS)
    if(BUILD_PERF_TESTS)
        add_subdirectory(PerfTests)
    endif()
    if (BUILD_PYTHON_BINDING)
        add_subdirectory(Python)
    endif()
endif()

#enabling to use HDIVPIOLA the contravariant Piola Mapping
option(USING_HDIVPIOLA "Whether the HDIVPIOLA is being used" OFF)
if(USING_HDIVPIOLA)
	#Adding a variable to hold this definition
	set(HDIVPIOLA HDIVPIOLA)
	add_definitions(-D${HDIVPIOLA})
endif(USING_HDIVPIOLA)

add_subdirectory(PerfUtil)
add_subdirectory(Util)
add_subdirectory(Common)
add_subdirectory(Save)
add_subdirectory(Integral)
add_subdirectory(LinearSolvers)
add_subdirectory(Matrix)
add_subdirectory(Topology)
add_subdirectory(Geom)
add_subdirectory(SpecialMaps)
add_subdirectory(Shape)
add_subdirectory(Refine)
add_subdirectory(External)
add_subdirectory(Material)
add_subdirectory(Mesh)
add_subdirectory(Analysis)
add_subdirectory(Multigrid)
add_subdirectory(Post)
add_subdirectory(Frontal)
add_subdirectory(StrMatrix)
add_subdirectory(Pre)
add_subdirectory(SubStruct)
add_subdirectory(Publications)
add_subdirectory("Random")
add_subdirectory(Optimization)
add_subdirectory(Exception)

#Creating groups for organized display of files into PZ modules (subfolders)
source_group(Common FILES ${pzcommon_STAT_SRCS})
source_group(Exception FILES ${pzexception_STAT_SRCS})
source_group(Save FILES ${pzsave_STAT_SRCS})
source_group(Util FILES ${pzutil_STAT_SRCS})
source_group(PerfUtil FILES ${pzperfutil_STAT_SRCS})
source_group(Matrix FILES ${pzmatrix_STAT_SRCS})
source_group("Linear Solvers" FILES ${pzlinearsolvers_STAT_SRCS})
source_group(Topology FILES ${pztopology_STAT_SRCS})
source_group(Geom FILES ${pzgeom_STAT_SRCS})
source_group("Special Maps" FILES ${pzspecialmaps_STAT_SRCS})
source_group(Mesh FILES ${pzmesh_STAT_SRCS})
source_group(Shape FILES ${pzshape_STAT_SRCS})
source_group(Refine FILES ${pzrefine_STAT_SRCS})
source_group(Material FILES ${pzmaterial_STAT_SRCS})
source_group(Analisys FILES ${pzanalysis_STAT_SRCS})
source_group(Multigrid FILES ${pzmultigrid_STAT_SRCS})
source_group(Post FILES ${pzpost_STAT_SRCS})
source_group(Frontal FILES ${pzfrontal_STAT_SRCS})
source_group("StructMatrix" FILES ${pzstrmatrix_STAT_SRCS})
source_group(Pre FILES ${pzpre_STAT_SRCS})
source_group("Substruct" FILES ${pzsubstruct_STAT_SRCS})
source_group(Integral FILES ${pzintegral_STAT_SRCS})
source_group(External FILES ${pzexternal_STAT_SRCS})
source_group(Publications FILES ${pzpublications_STAT_SRCS})
source_group(External\\Sloan FILES ${sloan_SRCS})
source_group(External\\FAD FILES ${FAD_SRCS})
source_group(External\\FAD\\Fad FILES ${FAD_fad_SRCS})
source_group(External\\FAD\\TinyFad FILES ${FAD_Tiny_SRCS})
source_group(External\\FAD\\TinyFadET FILES ${FAD_TinyET_SRCS})
source_group("Random" FILES ${pzrandom_STAT_SRCS})
source_group("Optimization" FILES ${pzoptimization_STAT_SRCS})

add_subdirectory(lib)


install(TARGETS pz 
		EXPORT PZTargets
		ARCHIVE
		DESTINATION ${NEOPZ_LIB})

if(BUILD_PROJECTS OR BUILD_TUTORIAL)
	add_subdirectory(Projects)
endif()

if(BUILD_UNITTESTING)
	add_subdirectory(UnitTest_PZ)
endif()

#Enabling tests
ENABLE_TESTING()
include(CTest)

#To insert cmakelist files as a group
if(BUILD_TUTORIAL)
	#Grouping CMakeList.txt for pz Projects - UnitTest_PZ
	file(GLOB_RECURSE makefiles *Lists.txt)

	#creating text files group
	source_group(CMakefilesPZ FILES ${makefilespz})
	source_group(CMakefilesProjects FILES ${makefilesproject})
	source_group(CMakefilesUnitTest FILES ${makefilesunittest})
	add_custom_target(CMAKE_Files SOURCES ${makefiles})
endif(BUILD_TUTORIAL)

#To insert doxygen documentation files
file(GLOB_RECURSE doxyfiles Documentation/*.h ./dox*.h)
add_custom_target(DOXYGEN_Files SOURCES  ${doxyfiles})

#To force XCode to use standard library c++ = libstdc++ (Required after upgrading to XCode 5)
#SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libstdc++")

set_target_properties(pz PROPERTIES XCODE_ATTRIBUTE_WARNING_CFLAGS "")


# Add all targets to the build-tree export set
export(TARGETS pz  FILE "${PROJECT_BINARY_DIR}/PZTargets.cmake")
 
# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE PZ)

set(INSTALL_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/CMake)

file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}"
   "${CMAKE_INSTALL_PREFIX}/${NEOPZ_INCLUDE}")   
   
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(PZConfig.cmake.in
  "${PROJECT_BINARY_DIR}/PZConfig.cmake" @ONLY)
  
# ... for the install tree
set(CONF_INCLUDE_DIRS "\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}"
	"/${REL_INCLUDE_DIR}/Python"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Util"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/PerfUtil"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Common"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Save"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Matrix"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Topology"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Geom"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/SpecialMaps"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Refine"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Shape"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Material"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Material/REAL"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Material/REAL/Plasticity"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Material/Complex"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Multigrid"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Mesh"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/StrMatrix"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Integral"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Frontal"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Pre"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Post"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Random"	
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Optimization"	
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Analysis"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/SubStruct"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/LinearSolvers"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/External"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/External/sloan"
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Publications"
)

if(USING_FAD)
	set(CONF_INCLUDE_DIRS ${CONF_INCLUDE_DIRS} 
            "\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/External/FAD"
            "\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/External/FAD/Fad"
            "\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/External/FAD/TinyFad"
            "\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/External/FAD/TinyFadET"
        )
endif(USING_FAD)

configure_file(PZConfig.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/PZConfig.cmake" @ONLY)
  
configure_file(PZConfig_Debug.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/PZConfig_Debug.cmake" COPYONLY)
   
# ... for both
configure_file(PZConfigVersion.cmake.in
  "${PROJECT_BINARY_DIR}/PZConfigVersion.cmake" @ONLY)
 
# Install the PZConfig.cmake and PZConfigVersion.cmake
install(FILES
        "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/PZConfig.cmake"
        "${PROJECT_BINARY_DIR}/PZConfigVersion.cmake"
        DESTINATION "${INSTALL_CMAKE_DIR}" )
  
install(FILES
        "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/PZConfig_Debug.cmake"
        DESTINATION "${INSTALL_CMAKE_DIR}" 
        CONFIGURATIONS Debug )  
 
# Install the export set for use with the install-tree
install(EXPORT PZTargets 
        DESTINATION "${INSTALL_CMAKE_DIR}")
