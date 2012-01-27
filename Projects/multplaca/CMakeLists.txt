
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(jornadas_SRCS jornadas.cpp)

add_executable(multplaca ${jornadas_SRCS})

target_link_libraries(multplaca ${KDE4_KDECORE_LIBS})

install(TARGETS multplaca ${INSTALL_TARGETS_DEFAULT_ARGS})

