########### next target ###############

file(GLOB header *.h ../TPZPlasticitySimulation.h)

file(GLOB source *.cpp ../TPZPlasticitySimulation.cpp)

add_executable(PlasticityCLI ${source} ${header})

target_link_libraries(PlasticityCLI pz)

