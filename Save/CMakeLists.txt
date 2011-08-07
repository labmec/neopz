# @file neopz/Save/CMakeLists.txt  -- CMake file for the Save module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)


set(pzsave_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


