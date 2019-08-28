
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS
    main.cpp
    pzadaptmesh.cpp
    pzcclonemesh.cpp
    pzgclonemesh.cpp)

add_executable(qualif ${main_SRCS})

target_link_libraries(qualif ${KDE4_KDECORE_LIBS})

install(TARGETS qualif ${INSTALL_TARGETS_DEFAULT_ARGS})

