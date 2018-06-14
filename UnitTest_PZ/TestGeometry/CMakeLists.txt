# @file neopz/UnitTest_PZ/TestGeometry/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	if (USING_FAD)

		include (CTestTestFile.cmake)

		set(pztestgeometry ${headers} ${sources} PARENT_SCOPE )

		add_executable(TestGeometry ${headers} ${sources})

		if(USING_BOOST)
			IF (WIN32)
				target_link_libraries(TestGeometry pz ${Boost_LIBRARIES})
			else()
				target_link_libraries(TestGeometry pz)
			endif(WIN32)
		else()
			target_link_libraries(TestGeometry pz)
		endif(USING_BOOST)

	endif(USING_FAD)

endif(BUILD_UNITTESTING)
