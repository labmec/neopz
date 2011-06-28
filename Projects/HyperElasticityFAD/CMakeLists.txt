########### next target ###############

set(main_SRCS main.cpp)

add_executable(HyperElasticity ${main_SRCS})

target_link_libraries(HyperElasticity pz pthread)

