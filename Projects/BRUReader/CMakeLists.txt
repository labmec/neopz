

set(main_SRCS main.cpp)

add_executable(BRUReader ${main_SRCS})

target_link_libraries(BRUReader pz)

install(TARGETS BRUReader ${INSTALL_TARGETS_DEFAULT_ARGS})


