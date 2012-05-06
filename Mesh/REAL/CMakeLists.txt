# @file neopz/Mesh/REAL/CMakeLists.txt  -- CMake file for REAL Materials

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzmesh_REAL_SRCS ${headers} ${sources} PARENT_SCOPE)