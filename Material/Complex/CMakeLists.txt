# @file neopz/Material/Complex/CMakeLists.txt  -- CMake file for Complex Materials

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzmaterial_COMPLEX_SRCS ${headers} ${sources} PARENT_SCOPE)