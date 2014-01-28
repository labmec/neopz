# @file neopz/PlaneFracture/CMakeLists.txt  -- CMake file for the SWXGraph module

file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(swxgraph_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)
