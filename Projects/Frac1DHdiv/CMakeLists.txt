file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(Frac1DHdiv ${main_SRCS})

target_link_libraries(Frac1DHdiv pz)

