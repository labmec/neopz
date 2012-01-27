
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS main.cpp)

add_executable(ShapeFAD ${main_SRCS})

target_link_libraries(ShapeFAD ${KDE4_KDECORE_LIBS} pz metis)

install(TARGETS ShapeFAD ${INSTALL_TARGETS_DEFAULT_ARGS})


