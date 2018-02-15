# @file neopz/UnitTest_PZ/TestPlasticity/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)

configure_file("ExpectedSandler1979.txt" "ExpectedSandler1979.txt" COPYONLY)

if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestplasticity ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestPlasticity ${headers} ${sources})

if(USING_BOOST)
	IF (WIN32)
		target_link_libraries(TestPlasticity pz ${Boost_LIBRARIES})
	ELSE()
		target_link_libraries(TestPlasticity pz)
	ENDIF()
else()
	target_link_libraries(TestPlasticity pz)
endif()

endif()

