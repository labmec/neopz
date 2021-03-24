file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(StabilizedHybridElasticity ${main_SRCS} ${headers})

target_link_libraries(StabilizedHybridElasticity pz)

