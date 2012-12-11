include_directories(..)

file(GLOB headers ../AnalyticalFunctions.h ../Elasticitytensordata.h ../ElasticMatInterface2D.h ../PoroElasticMatInterface2D.h ../MeshGeneration.h ../pzporoelastic2d.h ../tinyxml.h ../tinystr.h  )
file(GLOB sources ../AnalyticalFunctions.cpp ../Elasticitytensordata.cpp ../ElasticMatInterface2D.cpp ../PoroElasticMatInterface2D.cpp ../MeshGeneration.cpp ../pzporoelastic2d.cpp Validation1Dfinite.cpp Validation1DInfinite.cpp Validation2DKar.cpp Validation2DRad.cpp ValidationOilsubsidence.cpp ValidationMandelsProblem.cpp FaultReactivationTendency.cpp ValidationGenericCode.cpp ../tinystr.cpp ../tinyxml.cpp ../tinyxmlerror.cpp ../tinyxmlparser.cpp )

set(Mymain Main.cpp ${sources} ${headers})

add_executable(Poroelastic2D_Validation ${Mymain})
target_link_libraries(Poroelastic2D_Validation pz)

