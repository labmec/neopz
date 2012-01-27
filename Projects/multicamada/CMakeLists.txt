
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS multicamada.cpp)

add_executable(multicamada ${main_SRCS})

target_link_libraries(multicamada ${KDE4_KDECORE_LIBS} pz metis)

install(TARGETS multicamada ${INSTALL_TARGETS_DEFAULT_ARGS})


