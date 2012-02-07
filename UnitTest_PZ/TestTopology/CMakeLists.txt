# @file neopz/UnitTest_PZ/TestTopology/CMakeLists.txt  -- CMake file for unit test of the topology module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	set(pztesttopology ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestTopology ${headers} ${sources})

if(USING_BOOST)
	target_link_libraries(TestTopology pz boost_unit_test_framework)
else()
	target_link_libraries(TestTopology pz)
endif()

endif()

