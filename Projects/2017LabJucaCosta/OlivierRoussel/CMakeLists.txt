file(GLOB header *.h)
file(GLOB source *.cpp)

set(HPAdaptive_SRCS ${header} ${source})

add_executable(OlivierRoussel ${HPAdaptive_SRCS})

target_link_libraries(OlivierRoussel pz)

