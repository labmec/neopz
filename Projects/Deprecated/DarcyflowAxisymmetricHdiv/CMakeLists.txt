########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(DarcyflowAxiSymmetricHDiv ${sources} ${header})

target_link_libraries(DarcyflowAxiSymmetricHDiv pz)

