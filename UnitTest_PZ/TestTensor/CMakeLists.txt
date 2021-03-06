# @file neopz/UnitTest_PZ/TestMatrix/CMakeLists.txt  -- CMake file for unit test of tensor

if(USING_FAD AND BUILD_PLASTICITY_MATERIALS)
	if ((STATE_TYPE STREQUAL "double") AND (REAL_TYPE STREQUAL "double"))
		set(TestSources TestTensor.cpp)
		add_unit_test(TestTensor ${TestSources})
	endif()
endif()

