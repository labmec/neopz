# @file neopz/UnitTest_PZ/TestRando/CMakeLists.txt  -- CMake file for unit test of the random module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestrandom ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestRandom ${headers} ${sources})

if(USING_BOOST)
	IF (WIN32)
		target_link_libraries(TestRandom pz ${Boost_LIBRARIES})
	ELSE()
		target_link_libraries(TestRandom pz)
	ENDIF()
else()
	target_link_libraries(TestRandom pz)
endif()

endif()

