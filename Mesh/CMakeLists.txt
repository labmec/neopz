# @file neopz/Mesh/CMakeLists.txt  -- CMake file for the Mesh module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(pzmesh_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

