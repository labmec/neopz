
file(GLOB header *.h)

file(GLOB source *.cpp)

add_executable(HDivEstabilizado ${source} ${header})

target_link_libraries(HDivEstabilizado pz)

