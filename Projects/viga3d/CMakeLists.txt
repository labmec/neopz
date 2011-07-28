
########### next target ###############

set(main_SRCS viga3d.cpp)

add_executable(main ${main_SRCS})

target_link_libraries(main PZ)

install(TARGETS main ${INSTALL_TARGETS_DEFAULT_ARGS})


