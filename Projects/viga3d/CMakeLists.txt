
########### next target ###############

set(main_SRCS viga3d.cpp)

add_executable(viga3d ${main_SRCS})

target_link_libraries(viga3d pz pthread)

install(TARGETS viga3d ${INSTALL_TARGETS_DEFAULT_ARGS})


