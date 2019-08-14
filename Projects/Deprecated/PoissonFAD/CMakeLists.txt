
set(main_SRCS main.cpp)

add_executable(PoissonFAD ${main_SRCS})

target_link_libraries(PoissonFAD ${KDE4_KDECORE_LIBS} pz)

install(TARGETS PoissonFAD ${INSTALL_TARGETS_DEFAULT_ARGS})

