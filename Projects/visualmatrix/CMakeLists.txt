
########### next target ###############

set(main_SRCS main.cpp)

add_executable(visualmatrix ${main_SRCS})

target_link_libraries(visualmatrix pz pthread)

