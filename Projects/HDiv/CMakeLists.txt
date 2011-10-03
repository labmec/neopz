
file(GLOB header *.h)

set(main_SRCS HDiv.cpp)

add_executable(HDiv ${main_SRCS} ${header})

target_link_libraries(HDiv pz pthread)

