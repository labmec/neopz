
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(jornadas_SRCS jornadas.cpp)

add_executable(jornadas ${jornadas_SRCS})

target_link_libraries(jornadas ${KDE4_KDECORE_LIBS})

install(TARGETS jornadas ${INSTALL_TARGETS_DEFAULT_ARGS})

