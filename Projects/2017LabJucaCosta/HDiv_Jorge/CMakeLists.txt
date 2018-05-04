file(GLOB header *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(HDiv_Jorge ${main_SRCS} ${header})

target_link_libraries(HDiv_Jorge pz)

