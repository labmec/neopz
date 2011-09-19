
set(main_SRCS
    main.cpp
    pzadaptmesh.cpp
    pzcclonemesh.cpp
    pzgclonemesh.cpp)

add_executable(adaptivity ${main_SRCS})

target_link_libraries(adaptivity pz metis)

##install(TARGETS main ${INSTALL_TARGETS_DEFAULT_ARGS})

