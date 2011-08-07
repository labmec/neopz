# @file neopz/Cosys/CMakeLists.txt  -- CMake file for the Cosys module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(pzcosys_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

