file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(TrabalhosLNCC15 ${main_SRCS} ${headers})

target_link_libraries(TrabalhosLNCC15 pz)