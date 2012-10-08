# @file neopz/Projects/CursoPZ/Class_6/CMakeLists.txt  -- CMake file for the Class_5 of the CursoPZ tutorial - Elasticity Material

file(GLOB header *.h)

file(GLOB source *.cpp)

set(pzclass6 ${source} PARENT_SCOPE)

add_executable(CursoPZ_Elasticity2D ${source} ${header})

target_link_libraries(CursoPZ_Elasticity2D pz)

