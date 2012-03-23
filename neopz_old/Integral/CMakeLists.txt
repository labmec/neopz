# @file neopz/Integral/CMakeLists.txt  -- CMake file for the Integral module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzintegral_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)


