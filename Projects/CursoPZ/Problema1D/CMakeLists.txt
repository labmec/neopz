# @file neopz/Projects/CursoPZ/Problema1D/CMakeLists.txt  -- CMake file for the Class_5 of the CursoPZ tutorial - Elasticity Material

file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(Problema1D ${source} ${header})

target_link_libraries(Problema1D pz)

