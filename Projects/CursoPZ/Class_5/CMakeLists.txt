# @file neopz/Projects/CursoPZ/Class_5/CMakeLists.txt  -- CMake file for the Class_5 of the CursoPZ tutorial

file(GLOB source *.cpp)

set(pzclass5 ${source} PARENT_SCOPE)

add_executable(Class_5 ${source})

target_link_libraries(Class_5 pz pthread)

