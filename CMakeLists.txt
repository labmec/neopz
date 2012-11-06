# @file neopz/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 2.8.0)

project (PZ)

#Specify where the headers files must be finded
include_directories(
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

#To insert doxygen documentation file and other text files into the project environment (to facility changes)
#Grouping CMakeList.txt for pz Projects - UnitTest_PZ
file(GLOB_RECURSE makefiles *Lists.txt)

file(GLOB_RECURSE doxyfiles Documentation/*.h ./dox*.h)

#creating text files group 
source_group(CMakefilesPZ FILES ${makefilespz})
source_group(CMakefilesProjects FILES ${makefilesproject})
source_group(CMakefilesUnitTest FILES ${makefilesunittest})
add_custom_target(CMake_Files SOURCES ${makefiles})
add_custom_target(Doxygen_Files SOURCES  ${doxyfiles})

## inserting main definitions
#specifying the path to neopz source code
add_definitions(-DPZSOURCEDIR=\"${PROJECT_SOURCE_DIR}\")
#specify where the refinement patterns can be found (default behaviour when running the code)
add_definitions(-DREFPATTERNDIR=\"${PROJECT_SOURCE_DIR}/Refine/RefPatterns\")
#To check REAL as TPZFlopCount  - Jorge
#add_definitions(-Dcontar)

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
	# Edson: adding -DNDEBUG -DNODEBUG -DZERO_INTERNAL_RESIDU -DMAKEINTERNAL.
	# Philippe, why do we need ZERO_INTERNAL_RESIDU and MAKEINTERNAL?
 	#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -DNODEBUG -DZERO_INTERNAL_RESIDU -DMAKEINTERNAL" CACHE STRING "c++ release flags go here" FORCE)

#To search files or libraries
	set(SEARCH_DIRS ENV PATH /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/metis ./../externallibs/metis/include ./../externallibs/metis/lib ./..metis ./../metis/include ./../metis/lib ./externallibs/lib ./externallibs/include ./../externallibs/lib ./../externallibs/include)
#To install neopz/include
	set(NEOPZ_INCLUDE ${PROJECT_SOURCE_DIR}/../externallibs/neopz/include)
#To install neopz library
	set(NEOPZ_LIB ${PROJECT_SOURCE_DIR}/../externallibs/neopz/lib)

#ENABLING EXTERNAL LIBRARIES

#enabling pthread library
find_path(PTHREAD_INCLUDE NAMES PThread.h pthread.h PATHS ${SEARCH_DIRS})
find_library(PTHREAD_LIB NAMES libpthread.dylib pthread.lib pthreadVC2.lib libPThread.so libPThread.a libpthread.so libpthread.a PATHS ${SEARCH_DIRS})
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
			include_directories(${Boost_INCLUDE_DIRS})
			link_directories(${Boost_LIB_DIRS})
			link_libraries(${Boost_LIB_DIRS})
		ENDIF()
   	endif()
endif(USING_BOOST)

#enabling FAD library
option(USING_FAD "Whether the FAD package will be linked with the library" OFF)
if(USING_FAD)
	add_definitions(-D_AUTODIFF)
	find_path(FAD_INCLUDE fadtype.h fadType.h PATHS ${SEARCH_DIRS} ./../FAD ./FAD ./externallibs/FAD ./../externallibs/FAD)
	if(FAD_INCLUDE-NOTFOUND)
		set(FAD_INCLUDE CACHE PATH "Directory to the source of FAD")
	else()
		include_directories(${FAD_INCLUDE} ${FAD_INCLUDE}/Fad ${FAD_INCLUDE}/TinyFad ${FAD_INCLUDE}/TinyFadET)
	endif()
endif(USING_FAD)

#enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
if(USING_LOG4CXX)
	add_definitions(-DLOG4CXX)
	find_path(LOG4CXX_INCLUDE log4cxx PATHS ${SEARCH_DIRS} ./../log4cxx ./../log4cxx/include)
	find_library(LOG4CXX_LIB NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a PATHS ${SEARCH_DIRS} ./../log4cxx ./../log4cxx/lib ./../externallibs/lib/${CMAKE_CONFIGURATION_TYPE})
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
	find_path(TBB_INCLUDE tbb.h PATHS ${SEARCH_DIRS})
	find_library(TBB_LIB NAMES libtbb.dylib libtbb.so PATHS ${SEARCH_DIRS})
	if(TBB_INCLUDE-NOTFOUND)
		set (TBB_INCLUDE "" CACHE PATH "Directory where tbb header files can be found")
	else()
		include_directories(${TBB_INCLUDE})
	endif()
	if(TBB_LIB-NOTFOUND)
		set (TBB_LIB "" CACHE PATH "Directory where the tbb library can be found")
	else()
		link_libraries(${TBB_LIB})
	endif()
endif(USING_TBB)

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


#To building several examples for training and experimentation PZ
option(BUILD_TUTORIAL "To build examples for training or teaching PZ" OFF)
if(BUILD_TUTORIAL)
	add_definitions(-DBUILD_TUTORIAL)
endif()


#To building all the specific projects developed with PZ environment
option(BUILD_PROJECTS "To build specific projects." OFF)

#To building all the performance projects developed with PZ environment
option(BUILD_PERF_TESTS "To build the performance test." OFF)

#To select type of REAL
set(REALTypesAre "Types are: double, float, long double")
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
else()
  message (FATAL_ERROR "ERROR: Please specify a valid type for REAL. ${REALTypesAre}")
endif()
message ("Using REAL type = ${REAL_TYPE}")

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
message ("Using STATE type = ${STATE_TYPE}")

#Writing type definitions on config.h
if(WIN32)
	set(VC "#define VC")
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	#define use of pthread static lib.
	add_definitions(-DPTW32_STATIC_LIB)
endif()
set(CONFIG_H "${CMAKE_BINARY_DIR}/Common/config.h")
set(FILE_DEFINITIONS "${CONFIG_H}")
FILE(WRITE ${FILE_DEFINITIONS}
"#ifndef ${REAL_TYPE_DEF}\n#define ${REAL_TYPE_DEF}\n#endif\n"
"#ifndef ${STATE_TYPE_DEF}\n#define ${STATE_TYPE_DEF}\n#endif\n"
"${VC}\n"
)
#Specify where the config.h is located and add it to install list
include_directories("${CMAKE_BINARY_DIR}/Common")
install(FILES ${CONFIG_H} DESTINATION ${NEOPZ_INCLUDE})
set(pz_configh ${CONFIG_H})# PARENT_SCOPE)

#increment all the directories where the sources files is

add_subdirectory(PerfUtil)
if(BUILD_PERF_TESTS)
	add_subdirectory(PerfTests)
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


#################################################################3

#Creating log file to cmake process  -- it's temporary, only to check variable contents
#set(LogFile "${PROJECT_SOURCE_DIR}/LogCMake.txt"
#	CACHE STRING "Log file to control CMake process"
#	) 

#If your CMake version is less than 2.8.5 but up to 2.6 the following loop must to be commented, because the FIND command isn't in that 
#foreach(arq ${makefiles}) 
#	string(FIND ${arq} "Project" val)
#	if(${val} LESS 0)
#		string(FIND ${arq} "UnitTest" val)
#		if(${val} GREATER -1)
#			list(APPEND makefilesunittest ${arq})
#		else()
#			string(FIND ${arq} "lib" val)
#			if(val LESS 0)
#				string(FIND ${arq} "neopz/CMakeList" val)
#				if(val LESS 0)
#					list(APPEND makefilespz ${arq})
#				endif()
#			endif()
#		endif()
#	else()
#		list(APPEND makefilesproject ${arq})
#	endif()
#endforeach()
