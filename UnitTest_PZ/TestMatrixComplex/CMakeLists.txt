# @file neopz/UnitTest_PZ/TestMatrixComplex/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestmatrixcomplex ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestMatrixComplex ${headers} ${sources})

if(USING_BOOST)
	IF (WIN32)
		target_link_libraries(TestMatrixComplex pz ${Boost_LIBRARIES})
	ELSE()
		target_link_libraries(TestMatrixComplex pz)
	ENDIF()
else()
	target_link_libraries(TestMatrixComplex pz)
endif()

endif()

