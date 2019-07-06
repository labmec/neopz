
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS main.cc)

add_executable(LuisGuilherme ${main_SRCS})

target_link_libraries(LuisGuilherme ${KDE4_KDECORE_LIBS})

install(TARGETS LuisGuilherme ${INSTALL_TARGETS_DEFAULT_ARGS})

