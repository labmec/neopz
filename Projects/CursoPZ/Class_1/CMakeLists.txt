# @file neopz/Projects/CursoPZ/Class_1/CMakeLists.txt  -- CMake file for the Class_1: Numerical Integration of the CursoPZ tutorial

file(GLOB header *.h)

file(GLOB source *.cpp)

set(pzclass1 ${source} PARENT_SCOPE)

add_executable(CursoPZ1_NumIntegration ${source} ${header})

target_link_libraries(CursoPZ1_NumIntegration pz)

include (CTestTestFile.cmake)