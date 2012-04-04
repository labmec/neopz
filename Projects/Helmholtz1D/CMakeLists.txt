file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(Helmholtz1D ${main_SRCS})

target_link_libraries(Helmholtz1D pz pthread)
