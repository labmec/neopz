# @file neopz/SpecialMaps/CMakeLists.txt  -- CMake file for the SpecialMaps module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include/pz)

file(GLOB sources *.cpp)

set(pzspecialmaps_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)



