
file(GLOB header *.h)
file(GLOB source *.cpp)

add_executable(HDivEstabilizado ${source} ${headers})

target_link_libraries(HDivEstabilizado pz)

