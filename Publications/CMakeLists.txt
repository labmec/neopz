# @file neopz/Publications/CMakeLists.txt  -- CMake file for the Publications module


#ifndef BUILD_COMPLEX_PROJECTS ##projetos nao compativeis com state=complex<double>
file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzpublications_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)
#endif