########### next target ###############

add_executable(MHM main.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHM pz)

add_executable(MHMElast mainElast.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHMElast pz)

configure_file(porous.txt porous.txt COPYONLY)

add_executable(MHM_Simple MHMSimple.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHM_Simple pz)

