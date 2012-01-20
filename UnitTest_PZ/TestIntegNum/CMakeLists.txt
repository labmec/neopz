# @file neopz/UnitTest_PZ/TestMatrix/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	set(pztestintegnum ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestIntNum ${headers} ${sources})

	target_link_libraries(TestIntNum pz)

endif()

