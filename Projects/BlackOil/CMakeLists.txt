
set(main_SRCS
    main.cpp
    malha.cpp
    postprocess.cpp
    tpzarc3d.cpp
    tpzgeoblend.cpp)

add_executable(BlackOil ${main_SRCS})

target_link_libraries(main ${KDE4_KDECORE_LIBS})

install(TARGETS main ${INSTALL_TARGETS_DEFAULT_ARGS})

