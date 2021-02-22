# @file neopz/UnitTest_PZ/TestGeometry/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	if (USING_FAD)

		include (CTestTestFile.cmake)

		set(pztestgeometry ${headers} ${sources} PARENT_SCOPE )

		add_executable(TestGeometry ${headers} ${sources})
                if (WIN32)
                  target_link_libraries(TestGeometry pz ${Boost_LIBRARIES})
                else()
                  target_link_libraries(TestGeometry pz)
                  target_link_libraries(TestGeometry ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
                endif()

	endif(USING_FAD)

endif(BUILD_UNITTESTING)
