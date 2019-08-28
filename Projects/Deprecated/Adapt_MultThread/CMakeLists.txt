
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(Adapt_MultThread_SRCS
    main.cpp
    pzadaptmesh.cpp
    pzcclonemesh.cpp
    pzgclonemesh.cpp
    c0-simplequad.cpp
    c01-lshape-quad.cpp
    c02-lshape-triang.cpp
    c03-lshape-quadtri.cpp
    c04-simplehexa.cpp
    c05-cornercube.cpp
    c06-Tetra.cpp
    c07-Prism.cpp
    c08-Mixed.cpp
    c09-Plate-Cedric.cpp
    c10-Prism-Pref.cpp
    c11-Piram-Tetra.cpp
    c15-Exp-Hexa.cpp)

add_executable(Adapt_MultThread ${Adapt_MultThread_SRCS})

target_link_libraries(Adapt_MultThread ${KDE4_KDECORE_LIBS})

install(TARGETS Adapt_MultThread ${INSTALL_TARGETS_DEFAULT_ARGS})


