
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(HeMan_SRCS
    c0-simplequad.cpp
    c1-simplecube.cpp
    main.cc
    MeshReader.cpp
    Quad_1_Rib_Marked.cpp
    Quad_2_Rib_Ad_Marked.cpp)

add_executable(HeMan ${HeMan_SRCS})

target_link_libraries(HeMan ${KDE4_KDECORE_LIBS})

install(TARGETS HeMan ${INSTALL_TARGETS_DEFAULT_ARGS})

