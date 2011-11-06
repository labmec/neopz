# @file neopz/SpecialMaps/CMakeLists.txt  -- CMake file for the SpecialMaps module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzspecialmaps_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)



