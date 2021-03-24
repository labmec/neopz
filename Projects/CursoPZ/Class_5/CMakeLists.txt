# @file neopz/Projects/CursoPZ/Class_5/CMakeLists.txt  -- CMake file for the Class_5 of the CursoPZ tutorial

file(GLOB header *.h)

file(GLOB source *.cpp)

set(pzclass5 ${source} PARENT_SCOPE)

add_executable(CursoPZ5_ComputMesh ${source} ${header})

target_link_libraries(CursoPZ5_ComputMesh pz)

#include (CTestTestFile.cmake)
