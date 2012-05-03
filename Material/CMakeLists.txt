# @file neopz/Material/CMakeLists.txt  -- CMake file for the Material module

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

if(BUILD_COMPLEX_PROJECTS)
  add_subdirectory(Complex)
else()
  add_subdirectory(REAL)
endif(BUILD_COMPLEX_PROJECTS)

set(pzmaterial_STAT_SRCS ${headers} ${sources} ${pzmaterial_REAL_SRCS} ${pzmaterial_COMPLEX_SRCS} PARENT_SCOPE)