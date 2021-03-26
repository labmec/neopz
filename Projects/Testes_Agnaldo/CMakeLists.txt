file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(Testes_Agnaldo ${main_SRCS} ${headers})

target_link_libraries(Testes_Agnaldo pz)


