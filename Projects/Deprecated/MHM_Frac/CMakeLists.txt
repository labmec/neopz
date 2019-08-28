########### next target ###############

add_executable(FracFlow FracFlow.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(FracFlow pz)

add_executable(FracPre FracPre.cpp TPZFracSet.h TPZFracSet.cpp)

target_link_libraries(FracPre pz)

add_executable(FracFlowGmsh FracFlowGmsh.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp TPZFracSet.cpp TPZFracSet.h TPZFracSimulation.h TPZFracSimulation.cpp)

target_link_libraries(FracFlowGmsh pz)

