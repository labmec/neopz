# @file neopz/Projects/CursoPZ/Class_2/CMakeLists.txt  -- CMake file for the Class_2 of the CursoPZ tutorial

file(GLOB header *.h)

file(GLOB source *.cpp)

set(pzclass2 ${source} PARENT_SCOPE)

add_executable(CursoPZ_Class2 ${source} ${header})

target_link_libraries(CursoPZ_Class2 pz pthread)

