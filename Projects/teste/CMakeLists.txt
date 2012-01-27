
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS placa.cpp)

add_executable(test ${main_SRCS})

target_link_libraries(test ${KDE4_KDECORE_LIBS})

install(TARGETS test ${INSTALL_TARGETS_DEFAULT_ARGS})


	