# @file neopz/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 2.8.0)

project (PZ)

#Specify where the headers files must be finded
include_directories(
	"${PROJECT_SOURCE_DIR}/Python"
	"${PROJECT_SOURCE_DIR}/Util"
	"${PROJECT_SOURCE_DIR}/PerfUtil"
	"${PROJECT_SOURCE_DIR}/Common"
	"${PROJECT_SOURCE_DIR}/Save"
	"${PROJECT_SOURCE_DIR}/Matrix"
	"${PROJECT_SOURCE_DIR}/Topology"
	"${PROJECT_SOURCE_DIR}/Geom"
	"${PROJECT_SOURCE_DIR}/SpecialMaps"
	"${PROJECT_SOURCE_DIR}/Refine"
	"${PROJECT_SOURCE_DIR}/Shape"
	"${PROJECT_SOURCE_DIR}/Material"
	"${PROJECT_SOURCE_DIR}/Material/REAL"
	"${PROJECT_SOURCE_DIR}/Material/REAL/Plasticity"
	"${PROJECT_SOURCE_DIR}/Material/Complex"
	"${PROJECT_SOURCE_DIR}/Multigrid"
	"${PROJECT_SOURCE_DIR}/Mesh"
	"${PROJECT_SOURCE_DIR}/Mesh/REAL"
	"${PROJECT_SOURCE_DIR}/StrMatrix"
	"${PROJECT_SOURCE_DIR}/Integral"
	"${PROJECT_SOURCE_DIR}/Frontal"
	"${PROJECT_SOURCE_DIR}/Pre"
	"${PROJECT_SOURCE_DIR}/Post"
	"${PROJECT_SOURCE_DIR}/Analysis"
	"${PROJECT_SOURCE_DIR}/SubStruct"
	"${PROJECT_SOURCE_DIR}/LinearSolvers"
	"${PROJECT_SOURCE_DIR}/External" 
	"${PROJECT_SOURCE_DIR}/External/sloan"
)

## inserting main definitions
#specifying the path to neopz source code
set(PZSOURCEDIR ${PROJECT_SOURCE_DIR})
add_definitions(-DPZSOURCEDIR=\"${PZSOURCEDIR}\")
#specify where the refinement patterns can be found (default behaviour when running the code)
set(REFPATTERNDIR ${PROJECT_SOURCE_DIR}/Refine/RefPatterns)
add_definitions(-DREFPATTERNDIR=\"${REFPATTERNDIR}\")


#Defining variables
#To cmake
string(FIND ${CMAKE_CXX_FLAGS_DEBUG} "-DDEBUG" val)
if(${val} LESS 0)
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" CACHE STRING "c++ debug flags go here" FORCE)
endif()
string(FIND ${CMAKE_CXX_FLAGS_RELEASE} "-DNODEBUG" val)
if(${val} LESS 0)
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNODEBUG -DZERO_INTERNAL_RESIDU -DMAKEINTERNAL" CACHE STRING "c++ release flags go here" FORCE)
endif()

# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/pthread ./../externallibs/metis ./../externallibs/metis/include ./../externallibs/metis/lib ./..metis ./../metis/include ./../metis/lib ./externallibs/lib ./externallibs/include ./../externallibs/lib ./../externallibs/include /usr/lib/atlas-base/ /usr/include/atlas /usr/local/atlas/include /usr/local/atlas/lib )

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

#enabling Boost library 
option(USING_BOOST "Whether the Boost package will be linked with the library" OFF)
if(USING_BOOST)
	IF (WIN32)
		set(Boost_USE_STATIC_LIBS   ON)
		set(Boost_USE_MULTITHREADED ON)
		#add_definitions ( -DBOOST_ALL_NO_LIB )
	ENDIF()
	find_package( Boost COMPONENTS unit_test_framework REQUIRED )
	add_definitions(-DUSING_BOOST)
   	if(Boost_FOUND)
		IF (WIN32)
			include_directories(${Boost_INCLUDE_DIR})
			link_directories(${Boost_LIBRARY_DIRS})
			link_libraries(${Boost_LIBRARIES})
			#message ("WIN32 BoostINFO => IncludeDIR: ${Boost_INCLUDE_DIR}, LIBDIR: ${Boost_LIBRARY_DIRS}, LIBRARIES: ${Boost_LIBRARIES}")
		ELSE()
			include_directories(${Boost_INCLUDE_DIR})
			link_directories(${Boost_LIBRARY_DIRS})
			#link_libraries(${Boost_LIBRARIES})
			#message ("LINUX BoostINFO => IncludeDIR: ${Boost_INCLUDE_DIR}, LIBDIRs: ${Boost_LIBRARY_DIRS}, LIBRARIES: ${Boost_LIBRARIES}")
		ENDIF()
   	endif()
endif(USING_BOOST)


option(USING_FAD "Whether the FAD package will be linked with the library" OFF)
if(USING_FAD)
	#Adding a variable to hold this definition
	set(_AUTODIFF _AUTODIFF)
	add_definitions(-D${_AUTODIFF})
	include_directories(External/FAD External/FAD/Fad External/FAD/TinyFad External/FAD/TinyFadET)
endif(USING_FAD)

#enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
if(USING_LOG4CXX)
	#Adding a variable to hold this definition
	set(LOG4CXX LOG4CXX)
	add_definitions(-D${LOG4CXX})
	find_path(LOG4CXX_INCLUDE log4cxx PATHS ${SEARCH_DIRS} ${PROJECT_SOURCE_DIR}/../log4cxx ${PROJECT_SOURCE_DIR}/../log4cxx_src/include ${PROJECT_SOURCE_DIR}/../externallibs/log4cxx_src ${PROJECT_SOURCE_DIR}/../externallibs/log4cxx_src/include)
	find_library(LOG4CXX_LIB NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a PATHS ${SEARCH_DIRS} ./../log4cxx_src ./../log4cxx_src/lib ./../externallibs/log4cxx_src ./../externallibs/log4cxx_src/lib ./../externallibs/lib/${CMAKE_CONFIGURATION_TYPE})
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
	add_definitions(-DUSING_METIS)
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
	add_definitions(-DUSING_OPENSSL)
	link_libraries(crypto)
endif(USING_OPENSSL)

#enabling TBB library
option(USING_TBB "Whether the TBB library will be linked in" OFF)
if(USING_TBB)
	add_definitions(-DUSING_TBB)
	find_path(TBB_INCLUDE tbb.h PATHS ${SEARCH_DIRS}
	/opt/intel/tbb/include/tbb /usr/include/tbb /opt/local/include/tbb)
	find_library(TBB_LIB NAMES libtbb.dylib libtbb.so PATHS ${SEARCH_DIRS}  /opt/intel/tbb/lib )
	if(TBB_INCLUDE-NOTFOUND)
		set (TBB_INCLUDE "" CACHE PATH "Directory where tbb header files can be found")
	else()
		include_directories(/opt/local/include)
		include_directories(${TBB_INCLUDE})
	endif()
	if(TBB_LIB-NOTFOUND)
		set (TBB_LIB "" CACHE PATH "Directory where the tbb library can be found")
	else()
		link_libraries(${TBB_LIB})
	endif()
endif(USING_TBB)

#enabling LIKWID library
option(USING_LIKWID "Whether the likwid library will be linked in" OFF)
if(USING_LIKWID)
	add_definitions(-DUSING_LIKWID)
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
	add_definitions(-DUSING_LIBNUMA)
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

#enabling BLAS library
option(USING_BLAS "Whether the Basic Linear Algebra Subprograms (BLAS) library will be linked in" OFF)
if(USING_BLAS)
	add_definitions(-DUSING_BLAS)
	find_path(CBLAS_INCLUDE cblas.h PATHS ${SEARCH_DIRS})
	#find_library(ATLAS_LIB NAMES libatlas.dylib libtatlas.dylib atlas.lib libatlas.so libtatlas.so libatlas.a PATHS ${SEARCH_DIRS})
	find_library(CBLAS_LIB NAMES libcblas.dylib libtcblas.dylib cblas.lib libcblas.so libtcblas.so libcblas.a PATHS ${SEARCH_DIRS})
	if(CBLAS_INCLUDE-NOTFOUND)
		set (CBLAS_INCLUDE "" CACHE PATH "Directory where CBLAS header files can be found")
	else()
		include_directories(${CBLAS_INCLUDE})
	endif()
	#if(ATLAS_LIB-NOTFOUND)
	#	set (ATLAS_LIB "" CACHE PATH "Directory where the ATLAS library can be found")
	#endif()
	if(CBLAS_LIB-NOTFOUND)
		set (CBLAS_LIB "" CACHE PATH "Directory where the CBLAS library can be found")
	endif()
endif(USING_BLAS)

#enabling LAPACK library
option(USING_LAPACK "Whether the LAPACK library will be linked in" OFF)
if(USING_LAPACK)
	add_definitions(-DUSING_LAPACK)
	find_path(LAPACK_INCLUDE clapack.h PATHS ${SEARCH_DIRS})
	find_library(LAPACK_LIB NAMES libclapack.dylib PATHS ${SEARCH_DIRS})
	if(LAPACK_INCLUDE-NOTFOUND)
		set (LAPACK_INCLUDE "" CACHE PATH "Directory where LAPACK header files can be found")
	else()
		include_directories(${LAPACK_INCLUDE})
	endif()
	if(LAPACK_LIB-NOTFOUND)
		set (LAPACK_LIB "" CACHE PATH "Directory where the LAPACK library can be found")
	endif()
endif(USING_LAPACK)

#enabling the new skyline matrix format 
option(USING_NEW_SKYLMAT "Enable the new skyline matrix format" OFF)
if(USING_NEW_SKYLMAT)
        add_definitions(-DUSING_NEW_SKYLMAT)
endif(USING_NEW_SKYLMAT)

#enabling PAPI library 
option(USING_PAPI "Whether the PAPI library will be linked in" OFF)
if(USING_PAPI)
        add_definitions(-DUSING_PAPI)
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
        add_definitions(-DUSING_HWLOC)
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
        add_definitions(-DUSING_MKL)
        find_path(MKL_INCLUDE NAMES mkl.h PATHS ${SEARCH_DIRS} /opt/intel/mkl/include)
        find_library(MKL_LIB_INTEL NAMES libmkl_intel_lp64.so libmkl_intel_lp64.dylib PATHS
		${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/)
        find_library(MKL_LIB_CORE NAMES libmkl_core.so libmkl_core.dylib PATHS ${SEARCH_DIRS}
		/opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/)
        find_library(MKL_LIB_THREAD NAMES libmkl_intel_thread.so libmkl_intel_thread.dylib PATHS
		${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/)
find_library(COMPOSER_OMP NAMES libiomp5.so libiomp5.dylib PATHS ${SEARCH_DIRS}
/opt/intel/composer_xe/compiler/lib
/opt/intel/composerxe/lib/intel64
/opt/intel/lib)
#find_path(OMP_INCLUDE NAMES omp.h PATHS ${SEARCH_DIRS} /opt/intel/composer_xe_2013_sp1.1.103/compiler/include)
#include_directories(${OMP_INCLUDE})
        if(MKL_INCLUDE-NOTFOUND)
                set (MKL_INCLUDE "" CACHE PATH "Directory where mkl.h can be found")
        else()
                include_directories(${MKL_INCLUDE})
        endif()
        if(MKL_LIB-NOTFOUND)
                set (MKL_LIB_INTEL "" CACHE PATH "Directory where the mkl library can be found")
        endif()
endif(USING_MKL)


#To enable unit test using boost UnitTest framework
option(BUILD_UNITTESTING "Whether the Boost UnitTest will be linked in" OFF)
if(BUILD_UNITTESTING)
	add_definitions(-DBUILD_UNITTESTING)
	if(USING_BOOST)
		IF (NOT WIN32)
			link_libraries(${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
		ENDIF()
	endif()
endif()

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
	add_definitions(-DBUILD_TUTORIAL)
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
  add_definitions(-DREALdouble)
  set(REAL_TYPE_DEF "REALdouble")
elseif (REAL_TYPE STREQUAL "float")
  add_definitions(-DREALfloat)
  set(REAL_TYPE_DEF "REALfloat")
elseif (REAL_TYPE STREQUAL "long double")
  add_definitions(-DREALlongdouble)
  set(REAL_TYPE_DEF "REALlongdouble")
elseif (REAL_TYPE STREQUAL "pzfpcounter")
  add_definitions(-DREALpzfpcounter)
  set(REAL_TYPE_DEF "REALpzfpcounter")
else()
  message (FATAL_ERROR "ERROR: Please specify a valid type for REAL. ${REALTypesAre}")
endif()
message (STATUS "Using REAL type = ${REAL_TYPE}")

#To select type of STATE
set(STATETypesAre "Types are: double, float, long double, complex<float>, complex<double>, complex<long double>")
set(STATE_TYPE "double" CACHE STRING "${STATETypesAre}")

if (STATE_TYPE STREQUAL "double")
  add_definitions(-DSTATEdouble)
  set(STATE_TYPE_DEF "STATEdouble")
elseif (STATE_TYPE STREQUAL "float")
  add_definitions(-DSTATEfloat)
  set(STATE_TYPE_DEF "STATEfloat")
elseif (STATE_TYPE STREQUAL "long double")
  add_definitions(-DSTATElongdouble)
  set(STATE_TYPE_DEF "STATElongdouble")
elseif (STATE_TYPE STREQUAL "complex<float>")
  add_definitions(-DSTATEcomplexf)
  add_definitions(-DSTATE_COMPLEX)
  set(STATE_TYPE_DEF "STATEcomplexf")
  set (BUILD_COMPLEX_PROJECTS ON)
elseif (STATE_TYPE STREQUAL "complex<double>")
  add_definitions(-DSTATEcomplexd)
  add_definitions(-DSTATE_COMPLEX)
  set(STATE_TYPE_DEF "STATEcomplexd")
  set (BUILD_COMPLEX_PROJECTS ON)
elseif (STATE_TYPE STREQUAL "complex<long double>")
  add_definitions(-DSTATEcomplexld)
  add_definitions(-DSTATE_COMPLEX)
  set(STATE_TYPE_DEF "STATEcomplexld")
  set (BUILD_COMPLEX_PROJECTS ON)
else()
  message (FATAL_ERROR "ERROR: Please specify a valid type for STATE. ${STATETypesAre}")
endif()
message (STATUS "Using STATE type = ${STATE_TYPE}")

# Some definitions for windows
if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on config.h")
	add_definitions(-DVC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	#disabling All warnings
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W0")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W0")
	#define use of pthread static lib.
	add_definitions(-DPTW32_STATIC_LIB)
endif()

# Begin of config.h part
# This file contains all definitions used to compile the library pz
# this will be saved at include diretory in the install target
set (CONFIG_H "${PROJECT_BINARY_DIR}/Common/config.h")
configure_file ("${PROJECT_SOURCE_DIR}/Common/config.h.in" ${CONFIG_H})
# add it to include search path
include_directories("${PROJECT_BINARY_DIR}/Common")
# add the .h file in the install target
install (FILES ${CONFIG_H} DESTINATION ${NEOPZ_INCLUDE})
# End of config.h part

#increment all the directories where the sources files is
if(NOT BUILD_COMPLEX_PROJECTS)
    add_subdirectory(PerfUtil)
    if(BUILD_PERF_TESTS)
        add_subdirectory(PerfTests)
    endif()
    if (BUILD_PYTHON_BINDING)
        add_subdirectory(Python)
    endif()
endif()

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

add_subdirectory(lib)
if(BUILD_PROJECTS OR BUILD_TUTORIAL)
	add_subdirectory(Projects)
endif()
if(BUILD_UNITTESTING)
	add_subdirectory(UnitTest_PZ)
endif()

#Enabling tests
ENABLE_TESTING()

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
SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libstdc++")


