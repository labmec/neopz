file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain main.cpp ${sources} ${headers})

add_executable(Elasticity1D ${Mymain})

target_link_libraries(Elasticity1D pz)

