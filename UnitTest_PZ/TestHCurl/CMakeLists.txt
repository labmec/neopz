# @file neopz/UnitTest_PZ/TestBlend/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestmesh ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestHCurl ${headers} ${sources})
        if (WIN32)
	  target_link_libraries(TestHCurl pz ${Boost_LIBRARIES})
	else()
	  target_link_libraries(TestHCurl pz)
          target_link_libraries(TestHCurl ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
        endif()
endif(BUILD_UNITTESTING)
