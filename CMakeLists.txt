cmake_minimum_required (VERSION 2.8.5)

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
	"${PROJECT_SOURCE_DIR}/LinearSolvers"
	"${PROJECT_SOURCE_DIR}/External"
	"${PROJECT_SOURCE_DIR}/External/sloan"
	"${PROJECT_SOURCE_DIR}/PlaneFracture"
 )


#To insert doxygen documentation file and other text files into the project environment (to facility changes)
file(GLOB_RECURSE makefiles CMake*.txt README)
#Only CMakeList.txt corresponding any pz module but CMakeList.txt from UnitTest case can be considered
#If your CMake version is less than 2.8.5 but up to 2.6 the following loop must to be commented, because the FIND command isn't in that 
foreach(arq ${makefiles}) 
	string(FIND ${arq} "Project" val)
	if(${val} LESS 0)
		string(FIND ${arq} "UnitTest" val)
		if(${val} LESS 0)
			list(APPEND makefilescpy ${arq})
		endif()
	endif()
endforeach()

file(GLOB_RECURSE doxyfiles Documentation/*.h ./dox*.h)

#creating text files group 
source_group(CMakeFiles FILES ${makefilescpy})
source_group(DoxyFiles FILES ${doxyfiles})
add_custom_target(pzdocs SOURCES ${makefilescpy} ${doxyfiles})


#Defining variables
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" CACHE STRING "c++ debug flags go here" FORCE)

set (PTHREAD_INCLUDE "/usr/include" CACHE PATH "Directory where pthread.h can be found")
set (PTHREAD_LIB_INCLUDE "/usr/lib" CACHE PATH "Directory where the pthread library can be found")
include_directories(${PTHREAD_INCLUDE})
link_directories(${PTHREAD_LIB_INCLUDE})
link_libraries(pthread)

#specify where the refinement patterns can be found (default behaviour when running the code)
add_definitions(-DREFPATTERNDIR=\"${PROJECT_SOURCE_DIR}/Refine/RefPatterns\")

#enabling Boost library 
option(USING_BOOST "Whether the Boost package will be linked with the library" ON)

if(USING_BOOST)
	find_package( Boost )
   		if(Boost_FOUND)
			add_definitions(-DUSING_BOOST)
      		include_directories(${Boost_INCLUDE_DIRS})
			link_libraries(${Boost_LIB_DIRS})
   		endif()
endif()

#enabling FAD library
option(USING_FAD "Whether the FAD package will be linked with the library" OFF)

if(USING_FAD)
	find_path(FAD_INCLUDE FAD)
	if(NOT FAD_INCLUDE-NOTFOUND)
		set(FAD_INCLUDE CACHE PATH "Directory to the source of FAD")
		add_definitions(-D_AUTODIFF)
		include_directories(${FAD_INCLUDE} ${FAD_INCLUDE}/Fad ${FAD_INCLUDE}/TinyFad ${FAD_INCLUDE}/TinyFadET)
	endif()
endif(USING_FAD)

#enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
if(USING_LOG4CXX)
	set (LOG4CXX_INCLUDE "/usr/local/include" CACHE PATH "Directory where log4cxx header files can be found")
	set (LOG4CXX_LIB_INCLUDE "/usr/local/lib" CACHE PATH "Directory where the log4cxx library can be found")
	include_directories(${LOG4CXX_INCLUDE})
	link_directories(${LOG4CXX_LIB_INCLUDE})
	add_definitions(-DLOG4CXX -DPZSOURCEDIR=\"${PROJECT_SOURCE_DIR}\")
	link_libraries(log4cxx)
endif()

#To enable unit test using boost UnitTest framework
option(BUILD_UNITTESTING "Whether the Boost UnitTest will be linked in" OFF)
if(BUILD_UNITTESTING)
add_definitions(-DBUILD_UNITTESTING)
endif()
#To building several examples for training and experimentation PZ
option(BUILD_TUTORIAL "To build examples for training or teaching PZ" OFF)
if(BUILD_TUTORIAL)
add_definitions(-DBUILD_TUTORIAL)
endif()

if(BUILD_UNITTESTING)
	link_libraries(boost_unit_test_framework)
	add_subdirectory(UnitTest_PZ)
endif()

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
add_subdirectory(lib)

	add_subdirectory(PlaneFracture)
	add_subdirectory(Projects)

#Creating log file to cmake process  -- it's temporary, only to check variable contents
#set(LogFile "${PROJECT_SOURCE_DIR}/LogCMake.txt"
#	CACHE STRING "Log file to control CMake process"
#	) 
