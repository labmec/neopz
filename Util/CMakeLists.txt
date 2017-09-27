# @file neopz/Util/CMakeLists.txt  -- CMake file for the Util module

file(GLOB headers *.h Hash/*.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp Hash/*.cpp)

set(pzutil_STAT_SRCS ${headers} ${sources} ${pzutil_Hash_STAT_SRCS} PARENT_SCOPE)


