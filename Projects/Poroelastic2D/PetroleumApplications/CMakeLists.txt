include_directories(..)

file(GLOB headers ../AnalyticalFunctions.h  ../ElasticMatInterface2D.h ../MeshGeneration.h ../pzporoelastic2d.h ../tinyxml.h ../tinystr.h )
file(GLOB sources ../AnalyticalFunctions.cpp  ../ElasticMatInterface2D.cpp ../MeshGeneration.cpp ../pzporoelastic2d.cpp FaultReactivation.cpp ../tinystr.cpp ../tinyxml.cpp ../tinyxmlerror.cpp ../tinyxmlparser.cpp )

set(Mymain Main.cpp ${sources} ${headers})

add_executable(Poroelastic2D_App ${Mymain})
target_link_libraries(Poroelastic2D_App pz)