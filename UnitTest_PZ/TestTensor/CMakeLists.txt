# @file neopz/UnitTest_PZ/TestMatrix/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING AND USING_BOOST AND USING_FAD AND BUILD_PLASTICITY_MATERIALS)

	if ((STATE_TYPE STREQUAL "double") AND (REAL_TYPE STREQUAL "double"))
		include (CTestTestFile.cmake)

		set(pztesttensor ${headers} ${sources} PARENT_SCOPE )

		add_executable(TestTensor ${headers} ${sources})

		configure_file(tensor_and_eigensystem.txt tensor_and_eigensystem.txt COPYONLY)

		IF (WIN32)
			target_link_libraries(TestTensor pz ${Boost_LIBRARIES})
		ELSE()
		  target_link_libraries(TestTensor pz)
                  target_link_libraries(TestTensor ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
		ENDIF()
	endif()
endif()

