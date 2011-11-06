# @file neopz/Pre/CMakeLists.txt  -- CMake file for the Pre module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzpre_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

