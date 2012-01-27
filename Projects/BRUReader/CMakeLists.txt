
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS main.cpp)

add_executable(BRUReader ${main_SRCS})

target_link_libraries(BRUReader pz metis)

install(TARGETS BRUReader ${INSTALL_TARGETS_DEFAULT_ARGS})


