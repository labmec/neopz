
########### next target ###############

set(main_SRCS convhyper.cpp main.cpp)

add_executable(checkconv ${main_SRCS})

target_link_libraries(checkconv pz pthread)


