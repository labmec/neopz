########### next target ###############
file(GLOB sources *.cpp)

set(main_SRCS ${sources})

add_executable(Validacao ${main_SRCS})

target_link_libraries(Validacao pz pthread)

