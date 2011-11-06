# @file neopz/Shape/CMakeLists.txt  -- CMake file for the Shape module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzshape_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)




