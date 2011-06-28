########### next target ###############

set(main_SRCS main.cpp)

add_executable(matrix ${main_SRCS})

target_link_libraries(matrix pz pthread)

