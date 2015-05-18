########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(DarcyFlow3D ${sources} ${header})

target_link_libraries(DarcyFlow3D pz)

