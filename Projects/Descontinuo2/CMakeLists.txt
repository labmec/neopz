file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(Descontinuous2_Euler ${main_SRCS} ${headers})

target_link_libraries(Descontinuous2_Euler pz pthread)
