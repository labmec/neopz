########### next target ###############
file(GLOB header *.h)

set(MHMJorge_SRCS ${header} main.cpp TPZMatLaplacianHybrid.cpp meshgen.cpp HPAdaptiveProcesses.cpp)
set(MHMElastJorge_SRCS ${header} mainElast.cpp TPZMatLaplacianHybrid.cpp meshgen.cpp HPAdaptiveProcesses.cpp)

add_executable(MHM_Jorge ${MHMJorge_SRCS})

target_link_libraries(MHM_Jorge pz)

add_executable(MHMElast_Jorge ${MHMElastJorge_SRCS})

target_link_libraries(MHMElast_Jorge pz)

configure_file(porous.txt porous.txt COPYONLY)
