
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(hyperelastic_SRCS hyperelastic.cpp)

add_executable(materialhyperelastico ${hyperelastic_SRCS})

target_link_libraries(materialhyperelastico ${KDE4_KDECORE_LIBS})

install(TARGETS materialhyperelastico ${INSTALL_TARGETS_DEFAULT_ARGS})


