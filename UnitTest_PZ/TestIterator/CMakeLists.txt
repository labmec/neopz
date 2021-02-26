# @file neopz/UnitTest_PZ/TestIterator/CMakeLists.txt  -- CMake file for unit test of the TPZChunkVector iterator

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestiterator ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestIterator ${headers} ${sources})
        if (WIN32)
          target_link_libraries(TestIterator pz ${Boost_LIBRARIES})
        else()
          target_link_libraries(TestIterator pz)
          target_link_libraries(TestIterator ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
        endif()

endif()

