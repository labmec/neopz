# @file neopz/Topology/CMakeLists.txt  -- CMake file for the Topology module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pztopology_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)



