
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(RefPat_SRCS main.cpp)

add_executable(Test_RefPattern ${RefPat_SRCS})

target_link_libraries(Test_RefPattern ${KDE4_KDECORE_LIBS} pz metis)

install(TARGETS Test_RefPattern ${INSTALL_TARGETS_DEFAULT_ARGS})

