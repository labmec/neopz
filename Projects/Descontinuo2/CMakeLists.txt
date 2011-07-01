file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(EulerSimulation ${main_SRCS})

target_link_libraries(EulerSimulation pz pthread)
