file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(EulerDescontinuous2 ${main_SRCS} ${headers})

target_link_libraries(EulerDescontinuous2 pz pthread)
