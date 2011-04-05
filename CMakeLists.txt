cmake_minimum_required (VERSION 2.6)
project (PZ)

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
 )

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" CACHE STRING "c++ debug flags go here" FORCE)

set (PTHREAD_INCLUDE "/usr/include" CACHE PATH "directory where pthread.h can be found")
set (PTHREAD_LIB_INCLUDE "/usr/lib" CACHE PATH "directory where the pthread library can be found")
include_directories(${PTHREAD_INCLUDE})
link_directories(${PTHREAD_LIB_INCLUDE})
link_libraries(pthread)

option(USING_BOOST "Whether the Boost package will be linked with the library" ON)

if(USING_BOOST)
	
	find_package( Boost 1.36.0 )
   		if(Boost_FOUND)
			add_definitions(-DUSING_BOOST)
      		include_directories(${Boost_INCLUDE_DIRS})
   		endif()
endif()

option(USING_FAD "Whether the FAD package will be linked with the library" OFF)

if(USING_FAD)
	
	set(FAD_INCLUDE "/usr/local/include" CACHE PATH "directory to the source of FAD")
	add_definitions(-D_AUTODIFF)
    include_directories(${FAD_INCLUDE} ${FAD_INCLUDE}/Fad ${FAD_INCLUDE}/TinyFad ${FAD_INCLUDE}/TinyFadET)
endif()


option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
if(USING_LOG4CXX)
	set (LOG4CXX_INCLUDE "/usr/local/include" CACHE PATH "directory where log4cxx header files can be found")
	set (LOG4CXX_LIB_INCLUDE "/usr/local/lib" CACHE PATH "directory where the log4cxx library can be found")
	include_directories(${LOG4CXX_INCLUDE})
	link_directories(${LOG4CXX_LIB_INCLUDE})
	add_definitions(-DLOG4CXX -DPZSOURCEDIR=\"${PROJECT_SOURCE_DIR}\")
	link_libraries(log4cxx)
endif()

#configure_file(${PROJECT_SOURCE_DIR}/config.h.in ${PROJECT_BINARY_DIR}/config.h)
#configure_file(${PROJECT_SOURCE_DIR}/PZ_USER.mk.in ${PROJECT_BINARY_DIR}/PZ_USER.mk)

add_subdirectory(Util)
add_subdirectory(Common)
add_subdirectory(Save)
add_subdirectory(Integral)
#add_subdirectory(LinearSolvers)
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
#add_subdirectory(Cosys)
add_subdirectory(Multigrid)
add_subdirectory(Post)
add_subdirectory(Frontal)
add_subdirectory(StrMatrix)
add_subdirectory(Pre)
add_subdirectory(SubStruct)
add_subdirectory(lib)
add_subdirectory(Projects)

enable_testing()

add_test (TestShapeRuns Projects/testshape/TestShape)

add_test (TestGeomRuns Projects/testgeom/TestGeom)

include(CTest)
