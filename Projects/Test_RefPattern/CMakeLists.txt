
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(RefPat_SRCS main.cpp)

add_executable(RefPat ${RefPat_SRCS})

target_link_libraries(RefPat ${KDE4_KDECORE_LIBS} pz metis)

install(TARGETS RefPat ${INSTALL_TARGETS_DEFAULT_ARGS})

