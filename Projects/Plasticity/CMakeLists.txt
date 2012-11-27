file(GLOB headers *.h)
file(GLOB sources *.cpp)

file (GLOB to_remove_from_sources main.cpp mainGUI.cpp poroelastoplastic.cpp)
list (REMOVE_ITEM sources ${to_remove_from_sources})

set(main1_SRCS main.cpp poroelastoplastic.cpp ${sources} ${headers})
set(main2_SRCS mainGUI.cpp ${sources} ${headers})

if(USING_FAD)

	add_executable(Plasticity ${main1_SRCS})
	target_link_libraries(Plasticity pz)

	add_executable(PlasticityGUI ${main2_SRCS})
	target_link_libraries(PlasticityGUI pz)
	
	file(COPY SnubDodecahedron.txt DESTINATION .)

endif()

include (CTestTestFile.cmake)