# @file neopz/UnitTest_PZ/TestMatrix/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestmatrix ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestMatrix ${headers} ${sources})

if(USING_BOOST)
	IF (WIN32)
		target_link_libraries(TestMatrix pz ${Boost_LIBRARIES})
	ELSE()
		target_link_libraries(TestMatrix pz)
	ENDIF()
else()
	target_link_libraries(TestMatrix pz)
endif()

endif()

