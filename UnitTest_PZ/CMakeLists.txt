

if(BUILD_UNITTESTING)

	add_subdirectory(TestTopology)
	add_subdirectory(TestMatrix)
	add_subdirectory(TestMesh)

	source_group(UnitTestTopology FILES ${pztesttopology})
	source_group(UnitTestMatrix FILES ${pztestmatrix})
	source_group(UnitTestMesh FILES ${pztestmesh})
	
	set(pztest ${pztest} ${pztesttopology} ${pztestmatrix} ${pztestmesh} )

	add_custom_target(PZUnitTesting SOURCES ${pztest})


endif()

