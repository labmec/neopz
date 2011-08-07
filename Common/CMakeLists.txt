# @file neopz/Common/CMakeLists.txt  -- CMake file for the Common module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(pzcommon_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

