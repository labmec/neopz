# @file neopz/UnitTest_PZ/CMakeList.txt  -- First Cmake file for Unit Testing modules

if(BUILD_UNITTESTING)

set(RUNTIME_OUTPUT_DIRECTORY pz PARENT_SCOPE)

	add_subdirectory(TestTopology)
	add_subdirectory(TestMatrix)
	add_subdirectory(TestMesh)
	add_subdirectory(TestMaterial)

	source_group(UnitTestTopology FILES ${pztesttopology})
	source_group(UnitTestMatrix FILES ${pztestmatrix})
	source_group(UnitTestMesh FILES ${pztestmesh})
	
	#set(pztest ${pztest} ${pztesttopology} ${pztestmatrix} ${pztestmesh} )

	#add_custom_target(PZUnitTesting SOURCES ${pztest})

endif()

