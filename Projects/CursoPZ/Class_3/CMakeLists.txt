# @file neopz/Projects/CursoPZ/Class_3/CMakeLists.txt  -- CMake file for the Class_3 of the CursoPZ tutorial

file(GLOB header *.h)

file(GLOB source *.cpp)

set(pzclass3 ${source} PARENT_SCOPE)

add_executable(CursoPZ_Class3 ${source} ${header})

target_link_libraries(CursoPZ_Class3 pz pthread)

