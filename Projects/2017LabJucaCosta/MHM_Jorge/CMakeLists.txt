########### next target ###############
file(GLOB header *.h)
file(GLOB source *.cpp)

set(MHMJorge_SRCS ${header} ${source})

add_executable(MHM_Jorge ${MHMJorge_SRCS})

target_link_libraries(MHM_Jorge pz)

configure_file(porous_scaled.txt porous_scaled.txt COPYONLY)
