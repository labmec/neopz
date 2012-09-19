file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(LNCC1 ${main_SRCS} ${headers})

target_link_libraries(LNCC1 pz pthread)