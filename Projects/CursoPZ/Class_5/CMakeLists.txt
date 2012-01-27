# @file neopz/Projects/CursoPZ/Class_5/CMakeLists.txt  -- CMake file for the Class_5 of the CursoPZ tutorial

file(GLOB header *.h)

file(GLOB source *.cpp)

set(pzclass5 ${source} PARENT_SCOPE)

add_executable(CursoPZ_Class5 ${source} ${header})

target_link_libraries(CursoPZ_Class5 pz pthread)

