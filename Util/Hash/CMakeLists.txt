# @file neopz/Util/Hash/CMakeLists.txt  -- CMake file for the Hash module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Util/Hash)

file(GLOB sources *.cpp)

set(pzutil_Hash_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


