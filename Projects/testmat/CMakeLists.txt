
########### next target ###############

set(main_SRCS mainmat.cpp)

add_executable(testmat ${main_SRCS})

target_link_libraries(testmat pz pthread)

