# @file neopz/Projects/CursoPZ/Class_2/CMakeLists.txt  -- CMake file for the Class_2 of the CursoPZ tutorial

file(GLOB source *.cpp)

set(pzclass2 ${source} PARENT_SCOPE)

add_executable(Class_2 ${source})

target_link_libraries(Class_2 pz pthread)

