
file(GLOB header *.h)

set(main_SRCS HDiv_demo.cpp)

add_executable(HDiv_demo ${main_SRCS} ${header})

target_link_libraries(HDiv_demo pz)

