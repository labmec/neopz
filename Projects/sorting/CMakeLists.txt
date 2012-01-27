
include_directories(${KDE4_INCLUDES} ${KDE4_INCLUDE_DIR} ${QT_INCLUDES} )


########### next target ###############

set(sort_SRCS sort.cpp)

add_executable(sorting ${sort_SRCS})

target_link_libraries(sorting ${KDE4_KDECORE_LIBS} pz metis)

install(TARGETS sorting ${INSTALL_TARGETS_DEFAULT_ARGS})


