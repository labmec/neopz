# @file neopz/Mesh/CMakeLists.txt  -- CMake file for the Mesh module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
						   $<INSTALL_INTERFACE:.> 
						  )

file(GLOB headers *.h)

file(GLOB sources *.cpp)

if(BUILD_COMPLEX_PROJECTS)
  file (GLOB to_remove_from_headers pzhdivpressure.h pzflowcmesh.h tpzagglomeratemesh.h TPZAgglomerateEl.h TPZSBFemElementGroup.h TPZSBFemVolume.h)
  list (REMOVE_ITEM headers ${to_remove_from_headers})
  file (GLOB to_remove_from_sources pzflowcmesh.cpp TPZAgglomerateEl.cpp pzhdivpressure.cpp TPZSBFemElementGroup.cpp TPZSBFemVolume.cpp )
  list (REMOVE_ITEM sources ${to_remove_from_sources})
endif(BUILD_COMPLEX_PROJECTS)

if(NOT USING_MKL)
  file (GLOB to_remove_from_headers TPZSBFemElementGroup.h TPZSBFemVolume.h )
  list (REMOVE_ITEM headers ${to_remove_from_headers})
  file (GLOB to_remove_from_sources TPZSBFemElementGroup.cpp TPZSBFemVolume.cpp )
  list (REMOVE_ITEM sources ${to_remove_from_sources})
endif(NOT USING_MKL)

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Mesh)

set(pzmesh_STAT_SRCS ${headers} ${sources} ${pzmesh_REAL_SRCS} PARENT_SCOPE)

