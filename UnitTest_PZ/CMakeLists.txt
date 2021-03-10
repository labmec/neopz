# @file neopz/UnitTest_PZ/CMakeList.txt  -- First Cmake file for Unit Testing modules

if(BUILD_UNITTESTING)
	ENABLE_TESTING()

	set(RUNTIME_OUTPUT_DIRECTORY pz PARENT_SCOPE)

	if(USING_BOOST)

		add_subdirectory(TestFad)
		add_subdirectory(TestHCurl)
		add_subdirectory(TestTopology)
		add_subdirectory(TestMatrix)
		add_subdirectory(TestMesh)
		add_subdirectory(TestMHM)
		add_subdirectory(TestMaterial)
		add_subdirectory(TestIntegNum)
		add_subdirectory(TestIterator)
		add_subdirectory(TestGeometry)
		add_subdirectory(TestPlasticity)
		add_subdirectory(TestTensor)
		add_subdirectory(TestRandom)
		add_subdirectory(TestBlend)
		add_subdirectory(TestShape)
                add_subdirectory(TestMultithreading)

		source_group(UnitTestFad FILES ${pztestfad})
		source_group(UnitTestTopology FILES ${pztesttopology})
		source_group(UnitTestMatrix FILES ${pztestmatrix})
		source_group(UnitTestMesh FILES ${pztestmesh})
		source_group(UnitTestMHM FILES ${pztestmhm})
		source_group(UnitTestIntegNum FILES ${pztestintegnum})
		source_group(UnitTestIterator FILES ${pztestiterator})
		source_group(UnitTestGeometry FILES ${pztestgeometry})
		source_group(UnitTestPlasticity FILES ${pztestplasticity})
		source_group(UnitTestTensor FILES ${pztesttensor})
		source_group(UnitTestRandom FILES ${pztestrandom})
		source_group(UnitTestRandom FILES ${pztestshape})

	else()

		set(pztest ${pztest} ${pztestfad} ${pztesttopology} ${pztestmatrix} ${pztestmesh} ${pztestmhm} ${pztestintegnum}
				${pztestiterator} ${pztestrandom} ${pztestshape})

	endif()
endif()

