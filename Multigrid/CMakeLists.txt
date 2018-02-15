# @file neopz/Multigrid/CMakeLists.txt  -- CMake file for the multigrid module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Multigrid)

file(GLOB sources *.cpp)

set(pzmultigrid_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

