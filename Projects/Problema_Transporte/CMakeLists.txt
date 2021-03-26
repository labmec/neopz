
file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(Problema_Transporte ${source} ${header})

target_link_libraries(Problema_Transporte pz)

