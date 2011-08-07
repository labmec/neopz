# @file neopz/Util/CMakeLists.txt  -- CMake file for the Util module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(pzutil_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


