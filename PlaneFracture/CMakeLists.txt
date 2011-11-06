# @file neopz/PlaneFracture/CMakeLists.txt  -- CMake file for the PlaneFracture module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzplanefracture_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


