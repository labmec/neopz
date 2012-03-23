# @file neopz/Frontal/CMakeLists.txt  -- CMake file for the Frontal module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzfrontal_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

