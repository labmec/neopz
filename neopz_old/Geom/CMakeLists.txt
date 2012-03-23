# @file neopz/Geom/CMakeLists.txt  -- CMake file for the Geom module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzgeom_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

