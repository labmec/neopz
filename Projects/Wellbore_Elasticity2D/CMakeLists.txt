file(GLOB headers *.h)
file(GLOB sources *.cpp)

add_executable(Wellbore_Elasticity2D ${sources} ${headers})

target_link_libraries(Wellbore_Elasticity2D pz)

