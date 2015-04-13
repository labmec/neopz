# @file neopz/Publications/CMakeLists.txt  -- CMake file for the Publications module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzpublications_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)