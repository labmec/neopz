file(GLOB header *.h)
file(GLOB source *.cpp)

set(HPAdaptive_SRCS ${header} ${source})

add_executable(Adapt_HP_Jorge ${HPAdaptive_SRCS})

target_link_libraries(Adapt_HP_Jorge pz)

