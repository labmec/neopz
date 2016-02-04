#Project for Stokes problem

file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(Stokes ${main_SRCS})

target_link_libraries(Stokes pz)
