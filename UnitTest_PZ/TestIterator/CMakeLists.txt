# @file neopz/UnitTest_PZ/TestIterator/CMakeLists.txt  -- CMake file for unit test of the TPZChunkVector iterator

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestiterator ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestIterator ${headers} ${sources})

if(USING_BOOST)
	IF (WIN32)
		target_link_libraries(TestIterator pz ${Boost_LIBRARIES})
	ELSE()
		target_link_libraries(TestIterator pz)
	ENDIF()
else()
	target_link_libraries(TestIterator pz)
endif()

endif()

