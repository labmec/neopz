file(GLOB headers *.h)
file(GLOB sources *.cpp)

#REMOVER AS 3 SEGUINTES LINHAS QDO MOVER O MATERIAL PARA O PZ
include_directories(../MultiphysicsSimulation/)
set(mixedpoissonH ../MultiphysicsSimulation/mixedpoisson.h)
set(mixedpoissonCPP ../MultiphysicsSimulation/mixedpoisson.cpp)

set(main_SRCS ${sources} ${headers} ${mixedpoissonH} ${mixedpoissonCPP})

add_executable(ReducedSpace ${main_SRCS} ${headers})

target_link_libraries(ReducedSpace pz pthread)


