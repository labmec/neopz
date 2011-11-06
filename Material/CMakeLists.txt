# @file neopz/Material/CMakeLists.txt  -- CMake file for the Material module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzmaterial_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

