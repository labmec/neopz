########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(MHM ${sources} ${header})

target_link_libraries(MHM pz)

configure_file(porous.txt porous.txt COPYONLY)
