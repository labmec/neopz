
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS main.cpp)

add_executable(ReadTetGen ${main_SRCS})

target_link_libraries(ReadTetGen ${KDE4_KDECORE_LIBS})

install(TARGETS ReadTetGen ${INSTALL_TARGETS_DEFAULT_ARGS})


