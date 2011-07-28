
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(multigrid_SRCS main.cpp)

add_executable(multigrid ${multigrid_SRCS})

target_link_libraries(multigrid ${KDE4_KDECORE_LIBS} pthread)

install(TARGETS multigrid ${INSTALL_TARGETS_DEFAULT_ARGS})

