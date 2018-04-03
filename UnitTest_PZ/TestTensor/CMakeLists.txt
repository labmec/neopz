# @file neopz/UnitTest_PZ/TestMatrix/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztesttensor ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestTensor ${headers} ${sources})

        configure_file(tensor_and_eigensystem.txt tensor_and_eigensystem.txt COPYONLY)
if(USING_BOOST)
	IF (WIN32)
		target_link_libraries(TestTensor pz ${Boost_LIBRARIES})
	ELSE()
		target_link_libraries(TestTensor pz)
	ENDIF()
else()
	target_link_libraries(TestTensor pz)
endif()

endif()

