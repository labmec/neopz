# @file neopz/Post/CMakeLists.txt  -- CMake file for the Post module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzpost_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


