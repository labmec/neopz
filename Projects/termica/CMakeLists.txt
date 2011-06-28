########### next target ###############

set(main_SRCS termi.cpp fast.cpp)

add_executable(termica ${main_SRCS})

target_link_libraries(termica pz pthread)

