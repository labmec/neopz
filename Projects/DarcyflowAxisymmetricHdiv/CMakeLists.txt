########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(Darcyflow ${sources} ${header})

target_link_libraries(Darcyflow pz)

