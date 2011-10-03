
set(main_SRCS
    main.cc
    pzadaptmesh.cc
    pzcclonemesh.cc
    pzgclonemesh.cc)

add_executable(AdaptPhil ${main_SRCS})

target_link_libraries(AdaptPhil pz metis)

##install(TARGETS main ${INSTALL_TARGETS_DEFAULT_ARGS})

