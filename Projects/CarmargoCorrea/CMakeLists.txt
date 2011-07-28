
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(cccc_SRCS ctriangletests.cpp)

add_executable(cccc ${cccc_SRCS})

target_link_libraries(cccc ${KDE4_KDECORE_LIBS})

install(TARGETS cccc ${INSTALL_TARGETS_DEFAULT_ARGS})

