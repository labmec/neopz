
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(lara_SRCS
    MultiResMesh.cpp
    adapt.cpp
    main.cpp
    malhas.cpp
    pzexplfinvolanal.cpp
    pzexplfinvolanal.h
    tools.cpp)

add_executable(lara ${lara_SRCS})

target_link_libraries(lara ${KDE4_KDECORE_LIBS})

install(TARGETS lara ${INSTALL_TARGETS_DEFAULT_ARGS})


