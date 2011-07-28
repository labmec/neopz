
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(hyperelastic_SRCS hyperelastic.cpp)

add_executable(hyperelastic ${hyperelastic_SRCS})

target_link_libraries(hyperelastic ${KDE4_KDECORE_LIBS})

install(TARGETS hyperelastic ${INSTALL_TARGETS_DEFAULT_ARGS})


