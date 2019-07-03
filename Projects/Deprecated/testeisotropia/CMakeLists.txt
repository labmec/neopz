
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS testeisotropia.cpp)

add_executable(testeisotropia ${main_SRCS})

target_link_libraries(testeisotropia ${KDE4_KDECORE_LIBS})

install(TARGETS testeisotropia ${INSTALL_TARGETS_DEFAULT_ARGS})

