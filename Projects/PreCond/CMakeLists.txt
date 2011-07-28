
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(PreCond_SRCS PreCond.cpp)

add_executable(PreCond ${PreCond_SRCS})

target_link_libraries(PreCond ${KDE4_KDECORE_LIBS})

install(TARGETS PreCond ${INSTALL_TARGETS_DEFAULT_ARGS})

