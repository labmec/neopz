# @file neopz/UnitTest_PZ/TestFad/CMakeLists.txt  -- CMake file for unit test of FAD operations' consistency
file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING AND USING_BOOST AND USING_FAD AND BUILD_PLASTICITY_MATERIALS)

	if ((STATE_TYPE STREQUAL "double") AND (REAL_TYPE STREQUAL "double"))
		include (CTestTestFile.cmake)

		set(pztestfad ${headers} ${sources} PARENT_SCOPE )

		add_executable(TestFad ${headers} ${sources})

		if (WIN32)
			target_link_libraries(TestFad pz ${Boost_LIBRARIES})
		else()
		  target_link_libraries(TestFad pz)
		  target_link_libraries(TestFad ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
		endif()
	endif()
endif()
