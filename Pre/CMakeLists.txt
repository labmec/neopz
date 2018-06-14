# @file neopz/Pre/CMakeLists.txt  -- CMake file for the Pre module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)
file(GLOB sources *.cpp)

if (BUILD_COMPLEX_PROJECTS)
file (GLOB to_remove_from_headers pzreadmeshhr.h)
list (REMOVE_ITEM headers ${to_remove_from_headers})
file (GLOB to_remove_from_sources pzreadmeshhr.cpp)
list (REMOVE_ITEM sources ${to_remove_from_sources})
endif()

if(NOT USING_MKL OR BUILD_COMPLEX_PROJECTS)
  file (GLOB to_remove_from_headers TPZBuildSBFem.h )
  list (REMOVE_ITEM headers ${to_remove_from_headers})
  file (GLOB to_remove_from_sources TPZBuildSBFem.cpp )
  list (REMOVE_ITEM sources ${to_remove_from_sources})
endif(NOT USING_MKL OR BUILD_COMPLEX_PROJECTS)

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Pre)
set(pzpre_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)
