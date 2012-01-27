# @file neopz/Projects/CursoPZ/Class_4/CMakeLists.txt  -- CMake file for the Class_4 of the CursoPZ tutorial

file(GLOB header *.h)

file(GLOB source *.cpp)

set(pzclass4 ${source} PARENT_SCOPE)

add_executable(CursoPZ_Class4 ${source} ${header})

target_link_libraries(CursoPZ_Class4 pz pthread)

