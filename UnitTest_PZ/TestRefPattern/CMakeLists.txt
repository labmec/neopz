# @file neopz/UnitTest_PZ/TestBlend/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestgeometry ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestRefPattern ${headers} ${sources})

	if(USING_BOOST)
		IF (WIN32)
			target_link_libraries(TestRefPattern pz ${Boost_LIBRARIES})
		else()
			target_link_libraries(TestRefPattern pz)
		endif(WIN32)
	else()
		target_link_libraries(TestRefPattern pz)
	endif(USING_BOOST)

	

endif(BUILD_UNITTESTING)
