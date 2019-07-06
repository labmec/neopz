########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(LinearWave ${sources} ${header})

target_link_libraries(LinearWave pz)

