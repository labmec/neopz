# @file neopz/Matrix/CMakeLists.txt  -- CMake file for the Matrix module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzmatrix_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

