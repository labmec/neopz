# @file neopz/Util/CMakeLists.txt  -- CMake file for the Util module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzutil_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


