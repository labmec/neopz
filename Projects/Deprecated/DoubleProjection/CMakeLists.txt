file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(DoubleProjection ${main_SRCS} ${headers})

target_link_libraries(DoubleProjection pz)