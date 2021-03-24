########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(RefThread ${sources} ${header})

target_link_libraries(RefThread pz)

