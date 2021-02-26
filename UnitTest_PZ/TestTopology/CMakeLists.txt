# @file neopz/UnitTest_PZ/TestTopology/CMakeLists.txt  -- CMake file for unit test of the topology module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztesttopology ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestTopology ${headers} ${sources})

        if (WIN32)
	  target_link_libraries(TestTopology pz ${Boost_LIBRARIES})
	else()
	  target_link_libraries(TestTopology pz)
          target_link_libraries(TestTopology ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
        endif()

endif()

