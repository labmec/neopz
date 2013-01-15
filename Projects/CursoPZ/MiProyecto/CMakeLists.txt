# @file neopz/Projects/CursoPZ/MyProject/CMakeLists.txt  -- Simple project to test bugs in three dimensions

file(GLOB header *.h)

file(GLOB source *.cpp)

#set(pzclass6 ${source} PARENT_SCOPE)

add_executable(MiProyecto3D ${source} ${header})

target_link_libraries(MiProyecto3D pz)

