# @file neopz/CMakeList.txt  -- First file to read for CMake 

cmake_minimum_required (VERSION 2.8.0)

project (PZ)

#Specify where the headers files must be finded
include_directories(
	"${PROJECT_SOURCE_DIR}/Util"
	"${PROJECT_SOURCE_DIR}/Common"
	"${PROJECT_SOURCE_DIR}/Save"
	"${PROJECT_SOURCE_DIR}/Matrix"
	"${PROJECT_SOURCE_DIR}/Topology"
	"${PROJECT_SOURCE_DIR}/Geom"
	"${PROJECT_SOURCE_DIR}/SpecialMaps"
	"${PROJECT_SOURCE_DIR}/Refine"
	"${PROJECT_SOURCE_DIR}/Shape"
	"${PROJECT_SOURCE_DIR}/Material"
	"${PROJECT_SOURCE_DIR}/Multigrid"
	"${PROJECT_SOURCE_DIR}/Mesh"
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
	"${PROJECT_SOURCE_DIR}/PlaneFracture"
 )

#MESSAGE( STATUS "RUNTIME_DIR:         " ${RUNTIME_OUTPUT_DIRECTORY} )


#To insert doxygen documentation file and other text files into the project environment (to facility changes)
#Grouping CMakeList.txt for pz Projects - UnitTest_PZ
file(GLOB_RECURSE makefiles *Lists.txt)

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
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNODEBUG" CACHE STRING "c++ release flags go here" FORCE)
	endif()
#To search files or libraries
	set(SEARCH_DIRS ENV PATH /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/metis ./../externallibs/metis/include ./../externallibs/metis/lib ./..metis ./../metis/include ./../metis/lib)

#To install neopz/include
	set(NEOPZ_INCLUDE ${PROJECT_SOURCE_DIR}/../externallibs/neopz/include)
#To install neopz library
	set(NEOPZ_LIB ${PROJECT_SOURCE_DIR}/../externallibs/neopz/lib)

#ENABLING EXTERNAL LIBRARIES

#enabling pthread library
find_path(PTHREAD_INCLUDE NAMES PThread.h pthread.h PATHS ${SEARCH_DIRS})
find_path(PTHREAD_LIB NAMES libpthread.dylib pthread.lib libPThread.a libpthread.a PATHS ${SEARCH_DIRS})
if(PTHREAD_INCLUDE-NOTFOUND)
	set (PTHREAD_INCLUDE "/usr/include" CACHE PATH "Directory where pthread.h can be found")
endif()
if(PTHREAD_LIB-NOTFOUND)
	set (PTHREAD_LIB "/usr/lib" CACHE PATH "Directory where the pthread library can be found")
endif()
include_directories(${PTHREAD_INCLUDE})
link_directories(${PTHREAD_LIB})
link_libraries(pthread)

#enabling Boost library 
option(USING_BOOST "Whether the Boost package will be linked with the library" ON)
if(USING_BOOST)
	find_package( Boost )
   	if(Boost_FOUND)
		add_definitions(-DUSING_BOOST)
		include_directories(${Boost_INCLUDE_DIRS})
		link_directories(${Boost_LIB_DIRS})
		link_libraries(${Boost_LIB_DIRS})
   	endif()
endif(USING_BOOST)

#enabling FAD library
option(USING_FAD "Whether the FAD package will be linked with the library" OFF)
if(USING_FAD)
	find_path(FAD_INCLUDE fadtype.h PATHS ${SEARCH_DIRS} ./../FAD ./FAD)
	if(FAD_INCLUDE-NOTFOUND)
		set(FAD_INCLUDE CACHE PATH "Directory to the source of FAD")
	endif()
	add_definitions(-D_AUTODIFF)
	include_directories(${FAD_INCLUDE} ${FAD_INCLUDE}/Fad ${FAD_INCLUDE}/TinyFad ${FAD_INCLUDE}/TinyFadET)
endif(USING_FAD)

#enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" ON)
if(USING_LOG4CXX)
	find_path(LOG4CXX_INCLUDE log4cxx PATHS ${SEARCH_DIRS} ./../log4cxx ./../log4cxx/include)
	find_path(LOG4CXX_LIB NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.a liblog4cxx.so PATHS ${SEARCH_DIRS} ./../log4cxx ./../log4cxx/lib)
	if(LOG4CXX_INCLUDE-NOTFOUND)
		set (LOG4CXX_INCLUDE "/usr/local/include" CACHE PATH "Directory where log4cxx header files can be found")
	endif()
	if(LOG4CXX_LIB-NOTFOUND)
		set (LOG4CXX_LIB "/usr/local/lib" CACHE PATH "Directory where the log4cxx library can be found")
	endif()
	add_definitions(-DLOG4CXX)
	include_directories(${LOG4CXX_INCLUDE})
	link_directories(${LOG4CXX_LIB})
	link_libraries(log4cxx)
endif(USING_LOG4CXX)

#enabling METIS library
option(USING_METIS "Whether the METIS library will be linked in" OFF)
if(USING_METIS)
	find_path(METIS_INCLUDE metis.h PATHS ${SEARCH_DIRS})
	find_path(METIS_LIB NAMES libmetis.dylib metis.lib libmetis.a PATHS ${SEARCH_DIRS})
	if(METIS_INCLUDE-NOTFOUND)
		set (METIS_INCLUDE "/usr/local/include" CACHE PATH "Directory where metis header files can be found")
	endif()
	if(METIS_LIB-NOTFOUND)
		set (METIS_LIB "/usr/local/lib" CACHE PATH "Directory where the metis library can be found")
	endif()
	#verifying whether the metis is metis-5
	add_definitions(-DUSING_METIS)
	include_directories(${METIS_INCLUDE})
	link_directories(${METIS_LIB})
	link_libraries(metis)
endif(USING_METIS)

#To enable unit test using boost UnitTest framework
option(BUILD_UNITTESTING "Whether the Boost UnitTest will be linked in" OFF)
if(BUILD_UNITTESTING)
	add_definitions(-DBUILD_UNITTESTING)
	if(USING_BOOST)
		link_libraries(boost_unit_test_framework)
	endif()
endif()


#To building several examples for training and experimentation PZ
option(BUILD_TUTORIAL "To build examples for training or teaching PZ" OFF)
if(BUILD_TUTORIAL)
	add_definitions(-DBUILD_TUTORIAL)
endif()


#To building all the specific projects developed with PZ environment
option(BUILD_PROJECTS "To build specific projects." OFF)

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
  set(STATE_TYPE_DEF "STATEcomplexf")
  set (BUILD_COMPLEX_PROJECTS ON)
elseif (STATE_TYPE STREQUAL "complex<double>")
  add_definitions(-DSTATEcomplexd)
  set(STATE_TYPE_DEF "STATEcomplexd")
  set (BUILD_COMPLEX_PROJECTS ON)
elseif (STATE_TYPE STREQUAL "complex<long double>")
  add_definitions(-DSTATEcomplexld)
  set(STATE_TYPE_DEF "STATEcomplexld")
  set (BUILD_COMPLEX_PROJECTS ON)
else()
  message (FATAL_ERROR "ERROR: Please specify a valid type for STATE. ${STATETypesAre}")
endif()
message ("Using STATE type = ${STATE_TYPE}")

#Writing type definitions on config.h
set(FILE_DEFINITIONS "${CMAKE_SOURCE_DIR}/Common/config.h")
FILE(WRITE ${FILE_DEFINITIONS}
"#ifndef ${REAL_TYPE_DEF}\n#define ${REAL_TYPE_DEF}\n#endif\n"
"#ifndef ${STATE_TYPE_DEF}\n#define ${STATE_TYPE_DEF}\n#endif"
)

#increment all the directories where the sources files is
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
add_subdirectory(PlaneFracture)

add_subdirectory(lib)
add_subdirectory(Projects)
add_subdirectory(UnitTest_PZ)


#Creating log file to cmake process  -- it's temporary, only to check variable contents
#set(LogFile "${PROJECT_SOURCE_DIR}/LogCMake.txt"
#	CACHE STRING "Log file to control CMake process"
#	) 
