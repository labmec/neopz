file(GLOB header *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(HDiv ${main_SRCS} ${header})

target_link_libraries(HDiv pz)

