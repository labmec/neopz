# @file neopz/PerfUtil/CMakeLists.txt  -- CMake file for the PerfUtil module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzperfutil_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


