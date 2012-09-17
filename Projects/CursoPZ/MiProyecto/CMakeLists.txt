# @file neopz/Projects/CursoPZ/MyProject/CMakeLists.txt  -- CMake file for the Class_5 of the CursoPZ tutorial - Elasticity Material

file(GLOB header *.h)

file(GLOB source *.cpp)

#set(pzclass6 ${source} PARENT_SCOPE)

add_executable(MiProyecto2D ${source} ${header})

target_link_libraries(MiProyecto2D pz pthread)

