# @file neopz/Projects/CursoPZ/Class_3/CMakeLists.txt  -- CMake file for the Class_3 of the CursoPZ tutorial

file(GLOB source *.cpp)

set(pzclass3 ${source} PARENT_SCOPE)

add_executable(Class_3 ${source})

target_link_libraries(Class_3 pz pthread)

