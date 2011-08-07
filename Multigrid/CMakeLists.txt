# @file neopz/Multigrid/CMakeLists.txt  -- CMake file for the multigrid module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(pzmultigrid_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

