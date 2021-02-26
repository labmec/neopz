# @file neopz/UnitTest_PZ/TestPlasticity/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)

#configure_file("Sandler_Rubin_data_1979.txt" "Sandler_Rubin_data_1979.txt" COPYONLY)

if(BUILD_UNITTESTING AND USING_BOOST AND USING_FAD AND BUILD_PLASTICITY_MATERIALS)

	if ((STATE_TYPE STREQUAL "double") AND (REAL_TYPE STREQUAL "double"))
		include (CTestTestFile.cmake)

		set(pztestplasticity ${headers} ${sources} PARENT_SCOPE )

		add_executable(TestPlasticity ${headers} ${sources})

                if (WIN32)
                  target_link_libraries(TestPlasticity pz ${Boost_LIBRARIES})
                else()
                  target_link_libraries(TestPlasticity pz)
                  target_link_libraries(TestPlasticity ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
                endif()
	endif()
endif()

