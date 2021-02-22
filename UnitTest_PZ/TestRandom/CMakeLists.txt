# @file neopz/UnitTest_PZ/TestRando/CMakeLists.txt  -- CMake file for unit test of the random module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestrandom ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestRandom ${headers} ${sources})


	IF (WIN32)
	  target_link_libraries(TestRandom pz ${Boost_LIBRARIES})
	ELSE()
	  target_link_libraries(TestRandom pz)
          target_link_libraries(TestRandom ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
	ENDIF()

endif()

