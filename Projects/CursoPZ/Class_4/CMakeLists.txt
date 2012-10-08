# @file neopz/Projects/CursoPZ/Class_4/CMakeLists.txt  -- CMake file for the Class_4 of the CursoPZ tutorial

file(GLOB header *.h)

file(GLOB source *.cpp)

set(pzclass4 ${source} PARENT_SCOPE)

add_executable(CursoPZ4_GeoMesh ${source} ${header})

target_link_libraries(CursoPZ4_GeoMesh pz)

#include (CTestTestFile.cmake)
