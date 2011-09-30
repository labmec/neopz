file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(MultiphysicsSimulation ${main_SRCS})

target_link_libraries(MultiphysicsSimulation pz pthread)

