
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS
    main.cpp
    pzerror_ind.cpp
    tmbheidi.cpp
    tmbadapinterface.cpp
    tmbxubinshih.cpp)

add_executable(Error ${main_SRCS})

target_link_libraries(Error ${KDE4_KDECORE_LIBS})

install(TARGETS Error ${INSTALL_TARGETS_DEFAULT_ARGS})

