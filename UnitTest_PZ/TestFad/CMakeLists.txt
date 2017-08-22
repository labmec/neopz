# @file neopz/UnitTest_PZ/TestFad/CMakeLists.txt  -- CMake file for unit test of FAD operations' consistency
file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING AND USING_BOOST)

	include (CTestTestFile.cmake)

	set(pztestfad ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestFad ${headers} ${sources})

	if (WIN32)
		target_link_libraries(TestFad pz ${Boost_LIBRARIES})
	else()
		target_link_libraries(TestFad pz)
	endif()

endif()
