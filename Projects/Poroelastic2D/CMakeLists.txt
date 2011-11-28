file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(Poroelastic2D ${main_SRCS})

target_link_libraries(Poroelastic2D pz pthread)

