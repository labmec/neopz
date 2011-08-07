# @file neopz/Post/CMakeLists.txt  -- CMake file for the Post module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(pzpost_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


