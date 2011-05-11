
set(main_SRCS HDiv.cpp)

add_executable(HDiv ${main_SRCS})

target_link_libraries(HDiv pz pthread)

