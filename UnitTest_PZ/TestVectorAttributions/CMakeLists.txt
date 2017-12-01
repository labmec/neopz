# @file neopz/UnitTest_PZ/TestVectorAttributions/CMakeLists.txt  -- CMake file for unit test of the topology module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestvectorattributions ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestVectorAttributions ${headers} ${sources})

if(USING_BOOST)
	IF (WIN32)
		target_link_libraries(TestVectorAttributions pz ${Boost_LIBRARIES})
	ELSE()
		target_link_libraries(TestVectorAttributions pz)
	ENDIF()
else()
	target_link_libraries(TestVectorAttributions pz)
endif()

endif()

