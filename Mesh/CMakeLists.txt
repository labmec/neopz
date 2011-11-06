# @file neopz/Mesh/CMakeLists.txt  -- CMake file for the Mesh module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzmesh_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

