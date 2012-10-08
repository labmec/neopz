
set(main_SRCS main.cpp)

add_executable(AdptDescontinuo ${main_SRCS})

target_link_libraries(AdptDescontinuo ${KDE4_KDECORE_LIBS} pz)

install(TARGETS AdptDescontinuo ${INSTALL_TARGETS_DEFAULT_ARGS})

