# @file neopz/Mesh/CMakeLists.txt  -- CMake file for the Mesh module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

if(BUILD_COMPLEX_PROJECTS)
  file (GLOB to_remove_from_headers pzhdivpressure.h pzflowcmesh.h tpzagglomeratemesh.h TPZAgglomerateEl.h )
  list (REMOVE_ITEM headers ${to_remove_from_headers})
  file (GLOB to_remove_from_sources pzflowcmesh.cpp TPZAgglomerateEl.cpp pzhdivpressure.cpp )
  list (REMOVE_ITEM sources ${to_remove_from_sources})
endif(BUILD_COMPLEX_PROJECTS)

set(pzmesh_STAT_SRCS ${headers} ${sources} ${pzmesh_REAL_SRCS} PARENT_SCOPE)

