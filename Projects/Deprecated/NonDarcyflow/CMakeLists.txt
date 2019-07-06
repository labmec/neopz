########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(NonDarcyflow ${sources} ${header})

target_link_libraries(NonDarcyflow pz)

