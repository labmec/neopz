# @file neopz/Refine/CMakeLists.txt  -- CMake file for the Refine module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzrefine_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


add_subdirectory(RefPatterns)