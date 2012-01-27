# @file neopz/Projects/CursoPZ/Class_1/CMakeLists.txt  -- CMake file for the Class_1 of the CursoPZ tutorial

file(GLOB header *.h)

file(GLOB source *.cpp)

set(pzclass1 ${source} PARENT_SCOPE)

add_executable(CursoPZ_Class1 ${source} ${header})

target_link_libraries(CursoPZ_Class1 pz pthread)

