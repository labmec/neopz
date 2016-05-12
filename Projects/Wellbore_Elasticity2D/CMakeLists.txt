file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain main.cpp ${sources} ${headers})

add_executable(Wellbore_Elasticity2D ${Mymain})

target_link_libraries(Wellbore_Elasticity2D pz)

