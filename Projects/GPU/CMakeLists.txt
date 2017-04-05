#Project for GPU test

file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers})

add_executable(GPU ${main_SRCS})

target_link_libraries(GPU pz)


