# @file neopz/UnitTest_PZ/TestFad/CMakeLists.txt  -- CMake file for unit test of FAD operations' consistency

if (USING_FAD)
	set(TestSources FadUnitTest.cpp)
	add_unit_test(TestFad ${TestSources})
endif ()
