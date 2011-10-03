
file(GLOB headers *.h)

file(GLOB sources *.cpp)

add_executable(checkconv ${sources} ${headers})

target_link_libraries(checkconv pz pthread)


