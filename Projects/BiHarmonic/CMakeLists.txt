
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(main_SRCS main.cpp)

add_executable(BiHarmonic ${main_SRCS})

target_link_libraries(BiHarmonic ${KDE4_KDECORE_LIBS} pz metis pthread)

install(TARGETS BiHarmonic ${INSTALL_TARGETS_DEFAULT_ARGS})


