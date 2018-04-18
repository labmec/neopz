# @file neopz/Util/CMakeLists.txt  -- CMake file for the Util module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Util)

file(GLOB sources *.cpp)

add_subdirectory(Hash)

set(pzutil_STAT_SRCS ${headers} ${sources} ${pzutil_Hash_STAT_SRCS} PARENT_SCOPE)


