
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS main.cpp)

add_executable(NonLinBiharmonic ${main_SRCS})

target_link_libraries(NonLinBiharmonic ${KDE4_KDECORE_LIBS})

install(TARGETS NonLinBiharmonic ${INSTALL_TARGETS_DEFAULT_ARGS})

