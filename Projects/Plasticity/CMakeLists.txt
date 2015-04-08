file(GLOB headers *.h)
file(GLOB sources *.cpp)

file (GLOB to_remove_from_sources main.cpp mainGUI.cpp poroelastoplastic.cpp mainPV.cpp)
list (REMOVE_ITEM sources ${to_remove_from_sources})

set(main1_SRCS main.cpp poroelastoplastic.cpp ${sources} ${headers})
set(main3_SRCS mainPV.cpp)

#Plasticity works only with double type.
if ( STATE_TYPE STREQUAL "double" AND
     REAL_TYPE STREQUAL "double" AND 
     USING_FAD AND 
     BUILD_PLASTICITY_MATERIALS )

	  add_executable(Plasticity ${main1_SRCS})
	  target_link_libraries(Plasticity pz)

	  add_executable(PlasticityPV ${main3_SRCS})
	  target_link_libraries(PlasticityPV pz)

	#Note: you will need QT and QWT libraries to compile PlasticityQT
          #add_subdirectory(PlasticityQT)
          #add_subdirectory(interfaces_petro/WellboreQT)

	  file(COPY SnubDodecahedron.txt DESTINATION .)
	  file(COPY wellcil.txt DESTINATION .)
	  file(COPY ensaio_UCS_all_columns.txt DESTINATION .)
	  file(COPY ensaio_all_columns.txt DESTINATION .)

endif()

include (CTestTestFile.cmake)
