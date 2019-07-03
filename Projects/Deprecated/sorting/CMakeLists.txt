
set(sort_SRCS sort.cpp)

add_executable(sorting ${sort_SRCS})

target_link_libraries(sorting ${KDE4_KDECORE_LIBS} pz)

install(TARGETS sorting ${INSTALL_TARGETS_DEFAULT_ARGS})


