########### next target ###############

add_executable(MHM_Jorge main.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHM_Jorge pz)

add_executable(MHMElast_Jorge mainElast.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHMElast_Jorge pz)

configure_file(porous.txt porous.txt COPYONLY)
