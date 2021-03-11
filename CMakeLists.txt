# @file neopz/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.13.0)

#disable excessive warning messages by Microsoft VC compiler
#set (CMAKE_USER_MAKE_RULES_OVERRIDE "${CMAKE_CURRENT_LIST_DIR}/CompilerOptionsMSVC.cmake")

project(PZ LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#this allow us to define custom FindXXX.cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Definition of build type and REAL and STATE types
include(cmake/StandardPZSettings.cmake)

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


# if target_compile_options should be used in the following section,
# then maybe this could be useful:
# https://stackoverflow.com/questions/23995019/what-is-the-modern-method-for-setting-general-compile-flags-in-cmake/23995391
target_compile_definitions(pz PRIVATE
  "$<$<CONFIG:RELEASE>:NODEBUG ZERO_INTERNAL_RESIDU MAKEINTERNAL>"
  "$<$<CONFIG:DEBUG>:PZDEBUG>")

# Location to search files or libraries
set(EXTRA_SEARCH_DIRS ${EXTRA_SEARCH_DIRS})

#ENABLING EXTERNAL LIBRARIES
include(cmake/EnableThreads.cmake)
enable_threads(pz)

#enabling compiler interprocedural optimizations (IPO)
option(USING_IPO "Whether link and assembly support for compiler IPO will be activated or not (may cause linking to be very slow)." OFF)
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
	include(cmake/EnableBoost.cmake)
  enable_boost(pz)
endif()

option(USING_FAD "Whether the FAD package will be linked with the library" OFF)
if(USING_FAD)
	target_compile_definitions(pz PRIVATE _AUTODIFF)
endif()

#enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
if(USING_LOG4CXX)
  #this is the definition used in pz lib
  include(cmake/EnableLog4cxx.cmake)
  enable_log4cxx(pz)
  target_compile_definitions(pz PRIVATE LOG4CXX)
endif()

#enabling METIS library
option(USING_METIS "Whether the METIS library will be linked in" OFF)
if(USING_METIS)
    include(cmake/EnableMetis.cmake)
    enable_metis(pz)
endif()

#enabling OPENSSL library (used for MD5 computation)
option(USING_OPENSSL "Whether the OPENSSL library will be linked in" OFF)
if(USING_OPENSSL)
	target_link_libraries(pz PRIVATE crypto)
endif()

#enabling TBB library
option(USING_TBB "Whether the TBB library will be linked in" OFF)
if(USING_TBB)
    include(cmake/EnableTBB.cmake)
	enable_tbb(pz)
endif()

#enabling LIKWID library
option(USING_LIKWID "Whether the likwid library will be linked in" OFF)
if(USING_LIKWID)
    include(cmake/EnableLIKWID.cmake)
    enable_likwid(pz)
endif()

#enabling numa library
option(USING_LIBNUMA "Whether the Lib NUMA library will be linked in" OFF)
if(USING_LIBNUMA)
    include(cmake/EnableNUMA.cmake)
    enable_numa(pz)
endif()

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

#enabling MKL library
option(USING_MKL "Whether the MKL library will be linked in" OFF)
if(USING_MKL)
    include(cmake/EnableMKL.cmake)
    enable_mkl(pz)
endif(USING_MKL)

#enabling LAPACK library
option(USING_LAPACK "Whether the LAPACK library will be linked in" OFF)
if(USING_LAPACK)
    include(cmake/EnableLAPACK.cmake)
    enable_lapack(pz)
endif(USING_LAPACK)

#enabling BLAS library
option(USING_BLAS "Whether the Basic Linear Algebra Subprograms (BLAS) library will be linked in" OFF)
if(USING_BLAS)
    include(cmake/EnableBLAS.cmake)
    enable_blas(pz)
endif(USING_BLAS)

#enabling the new skyline matrix format
option(USING_NEW_SKYLMAT "Enable the new skyline matrix format" OFF)
if(USING_NEW_SKYLMAT)
		set(USING_NEW_SKYLMAT DUSING_NEW_SKYLMAT)
endif(USING_NEW_SKYLMAT)

#enabling PAPI library
option(USING_PAPI "Whether the PAPI library will be linked in" OFF)
if(USING_PAPI)
    include(cmake/EnablePAPI.cmake)
    enable_papi(pz)
endif(USING_PAPI)

#enabling Blaze
option(USING_BLAZE "Whether the Blaze package will be linked with the library" OFF)
if(USING_BLAZE)
    include(cmake/EnableBLAZE.cmake)
    enable_blaze(pz)
endif(USING_BLAZE)

#the following include is for making it easier to create pz targets
include(cmake/add_pz_target.cmake)

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
