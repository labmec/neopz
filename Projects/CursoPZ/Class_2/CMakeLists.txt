# @file neopz/Projects/CursoPZ/Class_2/CMakeLists.txt  -- CMake file for the Class_2 of the CursoPZ tutorial

file(GLOB header *.h)

file(GLOB source *.cpp)

set(pzclass2 ${source} PARENT_SCOPE)

add_executable(CursoPZ2_MatrixSolver ${source} ${header})

target_link_libraries(CursoPZ2_MatrixSolver pz)

include (CTestTestFile.cmake)

