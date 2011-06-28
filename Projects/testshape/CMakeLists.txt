
########### next target ###############

set(main_SRCS mainshape.cpp)

add_executable(TestShape ${main_SRCS})

target_link_libraries(TestShape pz pthread)

