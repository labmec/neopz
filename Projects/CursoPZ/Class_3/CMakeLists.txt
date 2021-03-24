# @file neopz/Projects/CursoPZ/Class_3/CMakeLists.txt  -- CMake file for the Class_3 of the CursoPZ tutorial

file(GLOB header *.h)

file(GLOB source *.cpp)

set(pzclass3 ${source} PARENT_SCOPE)

add_executable(CursoPZ3_Jacobian ${source} ${header})

target_link_libraries(CursoPZ3_Jacobian pz)

include (CTestTestFile.cmake)