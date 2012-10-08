
set(main_SRCS
    main.cc
    pzadaptmesh.cc
    pzcclonemesh.cc
    pzgclonemesh.cc)

add_executable(Adapt_Phil ${main_SRCS})

target_link_libraries(Adapt_Phil pz)
