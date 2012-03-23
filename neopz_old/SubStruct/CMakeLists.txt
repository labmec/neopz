# @file neopz/SubStruct/CMakeLists.txt  -- CMake file for the SubStruct module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzsubstruct_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


