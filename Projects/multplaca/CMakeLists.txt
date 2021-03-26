file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(multplaca ${sources} ${header})

target_link_libraries(multplaca pz)

#configure_file(multplaca.dat multplaca.dat COPYONLY)
#configure_file(multplaca2.dat multplaca2.dat COPYONLY)
#configure_file(multplaca3.dat multplaca3.dat COPYONLY)
#configure_file(onelayer.dat onelayer.dat COPYONLY)
