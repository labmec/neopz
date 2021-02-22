# @file neopz/UnitTest_PZ/TestIntegNum/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestintegnum ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestIntNum ${headers} ${sources})

        if (WIN32)
          target_link_libraries(TestIntNum pz ${Boost_LIBRARIES})
        else()
          target_link_libraries(TestIntNum pz)
          target_link_libraries(TestIntNum ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
        endif()
endif()

