# @file neopz/Material/REAL/Plasticity/CMakeLists.txt  -- CMake file for Plasticity Materials

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE})

file(GLOB sources *.cpp)

set(pzmaterial_Plasticity_SRCS ${headers} ${sources} PARENT_SCOPE)