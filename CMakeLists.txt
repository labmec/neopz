# @file neopz/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.13.0)

#disable excessive warning messages by Microsoft VC compiler
#set (CMAKE_USER_MAKE_RULES_OVERRIDE "${CMAKE_CURRENT_LIST_DIR}/cmake/CompilerOptionsMSVC.cmake")

project(PZ LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#this allow us to define custom FindXXX.cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Definition of build type and REAL and STATE types
include(cmake/StandardPZSettings.cmake)

## RPATH SETTINGS - BEGIN
# rpath settings must be set before any target is declared
include(GNUInstallDirs)

set(PZLIB_INSTALL_PREFIX "pz")
set(CMAKE_INSTALL_LIBDIR ${PZLIB_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_INSTALL_INCLUDEDIR ${PZLIB_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})
set(CMAKE_INSTALL_BINDIR ${PZLIB_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR})



file(RELATIVE_PATH relDir
${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)

if(APPLE)
set(basePoint @loader_path)
else()
set(basePoint $ORIGIN)
endif()
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH ${basePoint} ${basePoint}/${relDir})
##RPATH SETTINGS - END
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

# if target_compile_options should be used in the following section,
# then maybe this could be useful:
# https://stackoverflow.com/questions/23995019/what-is-the-modern-method-for-setting-general-compile-flags-in-cmake/23995391

#TODO: think if we don't need these definitions in RelWithDebugInfo builds
target_compile_definitions(pz PRIVATE
                           "$<$<CONFIG:RELEASE>:NODEBUG>"
                           "$<$<CONFIG:RELEASE>:ZERO_INTERNAL_RESIDU>"
                           "$<$<CONFIG:RELEASE>:MAKEINTERNAL>"
                           "$<$<CONFIG:DEBUG>:PZDEBUG>")
if(APPLE)                         
  set_target_properties(pz PROPERTIES XCODE_ATTRIBUTE_WARNING_CFLAGS "")
endif()
# Location to search files or libraries
set(EXTRA_SEARCH_DIRS ${EXTRA_SEARCH_DIRS})

#ENABLING EXTERNAL LIBRARIES

option(CMAKE_TBB_DEBUG "Outputs information when searching for the TBB library" OFF)
#this is for debugging any problems while finding mkl
option(CMAKE_MKL_DEBUG "Outputs information when searching for the MKL library" OFF)
mark_as_advanced(CMAKE_MKL_DEBUG)
mark_as_advanced(CMAKE_TBB_DEBUG)

include(cmake/EnableThreads.cmake)
enable_threads(pz)

IF(APPLE)
   # Mac OS X specific code
   find_library(ACCELERATE_LIB Accelerate)
   target_link_libraries(pz PRIVATE ${ACCELERATE_LIB})
ENDIF(APPLE)

#To enable unit test using Catch2 UnitTest framework
option(BUILD_UNITTESTING "Whether the Catch2 UnitTest will be linked in" OFF)
if(BUILD_UNITTESTING)
  include(cmake/add_unit_test.cmake)
  include(cmake/EnableCatch2.cmake)
  enable_catch2()
endif()


#enabling Boost library
option(USING_BOOST "Whether the Boost package will be linked with the library" OFF)
if(USING_BOOST)
	include(cmake/EnableBoost.cmake)
    enable_boost(pz)
endif()

#enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
if(USING_LOG4CXX)
  #this is the definition used in pz lib
  include(cmake/EnableLog4cxx.cmake)
  enable_log4cxx(pz)
  set(PZ_LOG ON) # Will add a define PZ_LOG on pz_config.h
  set(PZ_LOG4CXX_CONFIG_FILE ${PROJECT_SOURCE_DIR}/Util/log4cxx.cfg) #specify where log4cxx config is found
else()
  set(PZ_LOG OFF)
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
    target_compile_definitions(pz PRIVATE USING_OPENSSL)
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
include(cmake/add_pz_internal_target.cmake)

#To building several examples for training and experimentation PZ
option(BUILD_TUTORIAL "To build examples for training or teaching PZ" OFF)
if(BUILD_TUTORIAL)
	####add_definitions(-DBUILD_TUTORIAL)
endif()

#To building several examples for training and experimentation PZ
option(BUILD_BROKEN "To build projects which are broken" OFF)

#To building publication targets
option(BUILD_PUBLICATIONS "To build publications targets" OFF)
if(BUILD_PUBLICATIONS)
    add_subdirectory(Publications)
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
	set(VC ON) # Will add a define VC on pz_config.h
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	#disabling STL warning for std::_Vector_alloc when PZ_LOG is enabled
	#if (USING_PZ_LOG)
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
  find_package(Git QUIET)
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

# Provide a standard relative layout for the installation directory.
# Documentation directory is set to share/doc/PROJECT_NAME.


#this variable will contain all the public headers from neopz
set(PZ_PUBLIC_HEADERS "")
# add it to include search path
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
add_subdirectory(Random)
add_subdirectory(Optimization)
add_subdirectory(Exception)
#uncomment the following line to be scared
#message(${PZ_PUBLIC_HEADERS})
#############INSTALL#############

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  if(WIN32)
    set(CMAKE_INSTALL_PREFIX "C:/neopz" CACHE PATH "Directory where PZ will be installed" FORCE)
  else()
    set(CMAKE_INSTALL_PREFIX "/opt/neopz" CACHE PATH "Directory where PZ will be installed" FORCE)
  endif()
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# End of pz_config.h part

# Set a debug suffix automatically when PZ is built in Debug mode
set(CMAKE_DEBUG_POSTFIX "-debug")
set_target_properties(pz PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
install(TARGETS pz
        EXPORT NeoPZTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )
#where cmake config files are going to be set
set(PZ_CMAKE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/lib/cmake/neopz)

install(EXPORT NeoPZTargets
        FILE NeoPZTargets.cmake
        NAMESPACE NeoPZ::
        DESTINATION ${PZ_CMAKE_INSTALL_DIR}
        )
      
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/NeoPZConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/NeoPZConfig.cmake" @ONLY)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/NeoPZConfigVersion.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/NeoPZConfigVersion.cmake" @ONLY)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/NeoPZConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/NeoPZConfigVersion.cmake"
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/add_pz_target.cmake"
        DESTINATION ${PZ_CMAKE_INSTALL_DIR}
        )

include(cmake/InstallPZConfigFile.cmake)

#############EXTRA TARGETS#############
if(BUILD_PROJECTS OR BUILD_TUTORIAL)
    add_subdirectory(Projects)
endif()

if(BUILD_UNITTESTING)
    add_subdirectory(UnitTest_PZ)
    #Enabling tests
    ENABLE_TESTING()
    include(CTest)
endif()

if(NOT BUILD_COMPLEX_PROJECTS)
    if(BUILD_PERF_TESTS)
        add_subdirectory(PerfTests)
    endif()
endif()


#To build the documentation using doxygen
option(BUILD_DOCS "Whether to build the documentation using Doxygen" OFF)
#To build the documentation using doxygen+breathe+sphynx
option(BUILD_SPHINX_DOCS "Whether to build the documentation using Doxygen+Breathe+Sphinx" OFF)
if(BUILD_DOCS)
    add_subdirectory(docs)
elseif(BUILD_SPHINX_DOCS)
    message(FATAL_ERROR "For building sphinx docu BUILD_DOCS needs to be enabled")
endif()