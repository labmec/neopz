
set(main_SRCS main.cpp)

add_executable(ArtificialDiffusion ${main_SRCS})

target_link_libraries(ArtificialDiffusion pz metis)



