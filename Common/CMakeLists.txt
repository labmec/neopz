# @file neopz/Common/CMakeLists.txt  -- CMake file for the Common module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Common)

file(GLOB sources *.cpp)

set(pzcommon_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

