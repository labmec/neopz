# @file neopz/UnitTest_PZ/TestIntegNum/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestgeometry ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestGeometry ${headers} ${sources})

if(USING_BOOST)
	IF (WIN32)
		target_link_libraries(TestGeometry pz ${Boost_LIBRARIES})
	ELSE()
		target_link_libraries(TestGeometry pz)
	ENDIF()
else()
	target_link_libraries(TestGeometry pz)
endif()

endif()
