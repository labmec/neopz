
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(cccc_SRCS ctriangletests.cpp)

add_executable(CarmargoCorrea ${cccc_SRCS})

target_link_libraries(CarmargoCorrea ${KDE4_KDECORE_LIBS})

install(TARGETS CarmargoCorrea ${INSTALL_TARGETS_DEFAULT_ARGS})

