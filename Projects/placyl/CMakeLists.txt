
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS chapeu.c)

add_executable(placyl ${main_SRCS})

target_link_libraries(placyl ${KDE4_KDECORE_LIBS})

install(TARGETS placyl ${INSTALL_TARGETS_DEFAULT_ARGS})


