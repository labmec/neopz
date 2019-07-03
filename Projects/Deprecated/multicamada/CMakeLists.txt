
set(main_SRCS multicamada.cpp)

add_executable(multicamada ${main_SRCS})

target_link_libraries(multicamada ${KDE4_KDECORE_LIBS} pz)

install(TARGETS multicamada ${INSTALL_TARGETS_DEFAULT_ARGS})


