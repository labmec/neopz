# @file neopz/UnitTest_PZ/TestMatrix/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	set(pztestmatrix ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestMatrix ${headers} ${sources})

	target_link_libraries(TestMatrix pz boost_unit_test_framework)

endif()

