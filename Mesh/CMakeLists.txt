# @file neopz/Mesh/CMakeLists.txt  -- CMake file for the Mesh module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

if(BUILD_COMPLEX_PROJECTS)
else()
  add_subdirectory(REAL)
endif(BUILD_COMPLEX_PROJECTS)

set(pzmesh_STAT_SRCS ${headers} ${sources} ${pzmesh_REAL_SRCS} PARENT_SCOPE)

