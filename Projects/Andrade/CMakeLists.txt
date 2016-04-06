file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain Main.cpp ${sources} ${headers})

add_executable(Andrade ${Mymain})

target_link_libraries(Andrade pz)

