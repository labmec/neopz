# @file neopz/UnitTest_PZ/TestPlasticity/CMakeLists.txt  -- CMake file for plasticity unit test

if (USING_FAD AND BUILD_PLASTICITY_MATERIALS)
	if ((STATE_TYPE STREQUAL "double") AND (REAL_TYPE STREQUAL "double"))

		#configure_file("Sandler_Rubin_data_1979.txt" "Sandler_Rubin_data_1979.txt" COPYONLY)

		set(TestSources TestPlasticity.cpp)
		add_unit_test(TestPlasticity ${TestSources})

	endif ()
endif ()
