include_directories(..)

file(GLOB headers ../AnalyticalFunctions.h ../Elasticitytensordata.h ../ElasticMatInterface2D.h ../MeshGeneration.h ../pzporoelastic2d.h )
file(GLOB sources ../AnalyticalFunctions.cpp ../Elasticitytensordata.cpp ../ElasticMatInterface2D.cpp ../MeshGeneration.cpp ../pzporoelastic2d.cpp FaultReactivation.cpp)

set(Mymain Main.cpp ${sources} ${headers})

add_executable(Poroelastic2D_App ${Mymain})
target_link_libraries(Poroelastic2D_App pz pthread)