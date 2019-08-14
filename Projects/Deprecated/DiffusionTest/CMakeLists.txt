
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS main.cpp)

add_executable(DiffusionTest ${main_SRCS})

target_link_libraries(DiffusionTest ${KDE4_KDECORE_LIBS})

install(TARGETS DiffusionTest ${INSTALL_TARGETS_DEFAULT_ARGS})


