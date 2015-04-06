file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(TesteComplexo1D ${main_SRCS})

target_link_libraries(TesteComplexo1D pz)