
set(main_SRCS main.cpp)

add_executable(AltoDes ${main_SRCS})

target_link_libraries(AltoDes ${KDE4_KDECORE_LIBS} pz)

install(TARGETS AltoDes ${INSTALL_TARGETS_DEFAULT_ARGS})

