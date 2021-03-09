# @file neopz/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.13.0)

#disable excessive warning messages by Microsoft VC compiler
set (CMAKE_USER_MAKE_RULES_OVERRIDE "${CMAKE_CURRENT_LIST_DIR}/CompilerOptionsMSVC.cmake")

set(CMAKE_CONFIGURATION_TYPES Debug Release CACHE STRING  INTERNAL FORCE )

project(PZ LANGUAGES CXX)

#set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#this allow us to define custom FindXXX.cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
if (WIN32)
	if (NOT CMAKE_GENERATOR MATCHES "Ninja")
    set(CMAKE_INSTALL_PREFIX "C:/Program Files/PZ" CACHE PATH "Installation path, prepended onto install directories." FORCE) 
	endif()
else (WIN32)
  if ( NOT APPLE AND NOT CMAKE_BUILD_TYPE )
		set (CMAKE_BUILD_TYPE Release CACHE STRING "Options: Release Debug None MinSizeRel RelWithDebInfo" FORCE) 
	endif()
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


##TODOWIN32: should this change when on WIN32?
set(PZ_LIB_TYPE SHARED)


add_library(pz ${PZ_LIB_TYPE} "")

##TODOWIN32: check whether these distinctions are really necessary
##TODO: are they ok when private?
if (WIN32)
	# Prevents timespec redefinition problem with Visual Studio 2015
  target_compile_definitions(pz
    PRIVATE
    HAVE_STRUCT_TIMESPEC)
else (WIN32)
  target_compile_options(pz
    PRIVATE
    -fPIC
    -Wno-narrowing)
endif()

# Mac OS X specific code
if (APPLE)
  target_compile_definitions(pz PRIVATE MACOSX)
endif()

## inserting main definitions
#specifying the path to neopz source code
set(PZSOURCEDIR ${PROJECT_SOURCE_DIR})
#specify where the refinement patterns can be found (default behaviour when running the code)
set(REFPATTERNDIR ${PROJECT_SOURCE_DIR}/Refine/RefPatterns)

##TODO: need to refactor this following section
##maybe this could be useful https://stackoverflow.com/questions/23995019/what-is-the-modern-method-for-setting-general-compile-flags-in-cmake/23995391
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
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/local/opt/qt5/lib/cmake/ /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/metis ./../externallibs/metis/include ./../externallibs/metis/lib ./..metis ./../metis/include ./../metis/lib ./externallibs/lib ./externallibs/include ./../externallibs/lib ./../externallibs/include ./../externallibs/include/lapack /usr/lib/atlas-base/ /usr/include/atlas /usr/local/atlas/include /usr/local/atlas/lib )


#ENABLING EXTERNAL LIBRARIES

#prefer pthread library when there are other options available
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
#prefer -pthread compiler and linker flag
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
target_link_libraries(pz PRIVATE Threads::Threads)
#simple utility function for detecting if any linking is needed
#for using functions from <atomic>
include(cmake/check_for_atomic.cmake)
check_for_atomic(pz)
IF (USING_IPO)
  set_property(DIRECTORY PROPERTY INTERPROCEDURAL_OPTIMIZATION 1)
ENDIF()

IF(APPLE)
   # Mac OS X specific code
   find_library(ACCELERATE_LIB Accelerate)
   target_link_libraries(pz PRIVATE ${ACCELERATE_LIB})
ENDIF(APPLE)

#enabling Boost library
option(USING_BOOST "Whether the Boost package will be linked with the library" OFF)

#To enable unit test using boost UnitTest framework
option(BUILD_UNITTESTING "Whether the Boost UnitTest will be linked in" OFF)
if(BUILD_UNITTESTING)
	include(cmake/add_unit_test.cmake)
	set(USING_BOOST ON CACHE BOOL "Whether the Boost package will be linked with the library" FORCE)
endif()

if(USING_BOOST)
	find_package( Boost COMPONENTS unit_test_framework date_time atomic REQUIRED )
  target_link_libraries(pz PRIVATE Boost::date_time Boost::atomic)
endif(USING_BOOST)

option(USING_FAD "Whether the FAD package will be linked with the library" OFF)
if(USING_FAD)
	#Adding a variable to hold this definition
	set(_AUTODIFF _AUTODIFF)
endif(USING_FAD)

#enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
if(USING_LOG4CXX)
  #this is the definition used in pz lib
  set(LOG4CXX LOG4CXX)
	find_package(Log4cxx REQUIRED)
  target_link_libraries(pz PRIVATE ${Log4cxx_LIBRARY})
  target_include_directories(pz PRIVATE ${Log4cxx_INCLUDE_DIR})
  # if(WIN32) ##TODOWIN32: is it still needed?
	# target_link_libraries(pz PUBLIC odbc32.lib ws2_32.lib mswsock.lib)
  # endif()
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
  target_link_libraries(pz PRIVATE ${METIS_LIB})
  target_include_directories(pz PRIVATE ${METIS_INCLUDE})
endif(USING_METIS)

#enabling OPENSSL library (used for MD5 computation)
option(USING_OPENSSL "Whether the OPENSSL library will be linked in" OFF)
if(USING_OPENSSL)
	target_link_libraries(pz PRIVATE crypto)
endif(USING_OPENSSL)

#enabling TBB library
option(USING_TBB "Whether the TBB library will be linked in" OFF)
if(USING_TBB)
  find_package(TBB REQUIRED)
  if (NOT TBB_FOUND)
		message(FATAL_ERROR "Could not find TBB. If TBB is not needed, "
      "configure the project with -DUSING_TBB=OFF")
	endif()
  target_link_libraries(pz PRIVATE TBB::tbb)
  target_include_directories(pz PRIVATE ${TBB_INCLUDE_DIRS})
endif(USING_TBB)

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
		target_link_libraries(pz PRIVATE ${LIKWID_LIB})
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
		target_link_libraries(pz PRIVATE ${LIBNUMA_LIB})
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
  if(USING_MKL)
    set(BLA_VENDOR "Intel10_64lp") #tries to find MKL version
  endif()
  find_package(LAPACK REQUIRED)
  target_link_libraries(pz PRIVATE ${LAPACK_LIBRARIES} )
  foreach(LA_LIB LAPACK_LIBRARIES)
    if(${LA_LIB} MATCHES ".*mkl.*")
      message("Found LAPACK from MKL")
      target_compile_definitions(pz PRIVATE MKLLAPACK)
      #USING_BLAS=ON if USING_LAPACK=ON in pz_config.h.in
      target_compile_definitions(pz PRIVATE MKLBLAS)
    endif()
  endforeach()
endif(USING_LAPACK)

#enabling BLAS library
if(USING_BLAS)
  if(USING_MKL)
    set(BLA_VENDOR "Intel10_64lp") #tries to find MKL version
  endif()
  # set(BLA_VENDOR "GENERIC")
  find_package(BLAS REQUIRED)
  #target_link_libraries(pz BLAS::BLAS) this is ok in cmake 3.18 onwards
  target_link_libraries(pz PRIVATE ${BLAS_LIBRARIES})
  foreach(BL_LIB BLAS_LIBRARIES)
    if(${BL_LIB} MATCHES ".*mkl.*")
      message("Found BLAS from MKL")
      target_compile_definitions(pz PRIVATE MKLBLAS)
    endif()
  endforeach()
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
  target_link_libraries(pz PRIVATE ${PAPI_LIB})
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
  target_link_libraries(pz PRIVATE ${HWLOC_LIB})
endif(USING_HWLOC)



#enabling MKL library
option(USING_MKL "Whether the MKL library will be linked in" OFF)
if(USING_MKL)
  target_compile_definitions(pz PRIVATE MKLBLAS)
  target_compile_definitions(pz PRIVATE MKLLAPACK)
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
    target_link_libraries(pz PRIVATE ${MKL_LIB_INTEL})
    target_link_libraries(pz PRIVATE ${MKL_LIB_CORE})
    target_link_libraries(pz PRIVATE ${MKL_LIB_THREAD})
    target_link_libraries(pz PRIVATE ${COMPOSER_OMP})
  endif()
  target_link_libraries(pz PRIVATE ${MKL_LIB_INTEL} ${MKL_LIB_CORE} ${MKL_LIB_THREAD} ${COMPOSER_OMP})
  target_include_directories(pz PRIVATE ${MKL_INCLUDE})  
endif(USING_MKL)

#enabling Blaze
option(USING_BLAZE "Whether the Blaze package will be linked with the library" OFF)
if(USING_BLAZE)
        if ( NOT USING_MKL AND NOT USING_LAPACK )
                message( FATAL_ERROR "Blaze library requires either USING_MKL or USING_LAPACK enabled !")
        endif() 
	find_package( blaze REQUIRED)
	if( blaze_FOUND )
		#link_libraries(${MKL_LIB_INTEL})
	    add_library( blaze_target INTERFACE )
    else()
        set (USING_BLAZE "" CACHE PATH "Directory where blaze can be found")
    endif()
    target_compile_definitions(pz PRIVATE USING_BLAZE)
    target_link_libraries(pz PRIVATE blaze::blaze) # Blaze docs enforce this line. Mandatory under Windows. In future, should include "INTERFACE" modifier.
endif(USING_BLAZE)

#the following include is for making it easier to create pz targets
include(cmake/add_pz_target.cmake)


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
	add_definitions(-DPTW32_STATIC_LIB)
        #enabling /bigobj
        add_definitions("/bigobj")
        if (NOT CMAKE_GENERATOR MATCHES "Ninja")
            add_definitions("/MP")
        endif()
        if (CMAKE_GENERATOR MATCHES "NMake")
            set(CMAKE_MAKE_PROGRAM D:/Downloads/jom_1_1_3/jom.exe CACHE STRING "Path to jom.exe" FORCE)
            message("CMAKE_MAKE_PROGRAM set to 'D:/Downloads/jom_1_1_3/jom.exe' to enable parallel builds in Windows")
            FOREACH(FLAG_TYPE EXE MODULE SHARED)
                # (NMake only) No incremental linking; .pdb file must be generated or QtCreator is not able to debug (via CDB debugger)
                SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "/DEBUG /INCREMENTAL:NO" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
                MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
            ENDFOREACH ()
        else()
            FOREACH(FLAG_TYPE EXE MODULE SHARED)
                # (MSVC and Ninja only)No incremental linking; preventing .pdb unnecessary file generation overhead: no /DEBUG flag in linker flags
                SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "/INCREMENTAL:NO" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
                MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
            ENDFOREACH ()
        endif()
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

# INSTALL part
set(NEOPZ_INSTALL_DIR pzlib)
#To install neopz/include
set(NEOPZ_INCLUDE ${NEOPZ_INSTALL_DIR}/include)
#To install neopz library
set(NEOPZ_LIB ${NEOPZ_INSTALL_DIR}/lib)

install (FILES ${CONFIG_H} DESTINATION ${NEOPZ_INCLUDE} COMPONENT libraries)
# End of pz_config.h part

#increment all the directories where the sources files is
if(NOT BUILD_COMPLEX_PROJECTS)
    if(BUILD_PERF_TESTS)
        add_subdirectory(PerfTests)
    endif()
    if (BUILD_PYTHON_BINDING)
        add_subdirectory(Python)
    endif()
endif()

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

#TODOWIN32: Visual Studio source groups (filters) used to be generated for every subdirectory as follows:
# source_group(Common FILES ${pzcommon_STAT_SRCS})
# Currently, we don't have these source variables. If the source groups are really needed in Visual Studio,
# we need to think of an alternative way of creating them.

install(TARGETS pz 
		EXPORT PZTargets
		LIBRARY DESTINATION ${NEOPZ_LIB})

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

set(INSTALL_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/${NEOPZ_INSTALL_DIR}/CMake)

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
	"\${PZ_CMAKE_DIR}/${REL_INCLUDE_DIR}/Exception"
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

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
