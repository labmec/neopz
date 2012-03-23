# @file neopz/Analysis/CMakeLists.txt  -- CMake file for Analysis module. 

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzanalysis_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

