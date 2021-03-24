########### next target ###############

file(GLOB header *.h)

file(GLOB sources *.cpp)

add_executable(Gridgeneration ${sources} ${header})

target_link_libraries(Gridgeneration pz) 

