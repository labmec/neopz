file(GLOB header *.h)
file(GLOB source *.cpp)

set(HPAdaptive_SRCS ${header} ${source})

add_executable(HPAdaptive ${HPAdaptive_SRCS})

target_link_libraries(HPAdaptive pz)

