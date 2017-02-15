
file(GLOB header *.h)

set(main_SRCS HDiv.cpp)

add_executable(HDivGirkmann ${main_SRCS} ${header})

target_link_libraries(HDivGirkmann pz)

