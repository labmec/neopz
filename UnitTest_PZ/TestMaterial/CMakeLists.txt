# @file neopz/UnitTest_PZ/TestMaterial/CMakeLists.txt  -- CMake file for unit test of the matrix module

file(GLOB headers *.h)
file(GLOB sources *.cpp)

#file(COPY "CubeStiffMatrix.txt" DESTINATION ../../../NeoPZ_CMake/UnitTest_PZ/TestMaterial/Debug)
file(COPY "CubeStiffMatrix.txt" DESTINATION Debug)

if(BUILD_UNITTESTING)

	set(pztestmaterial ${headers} ${sources} PARENT_SCOPE )

	add_executable(TestMaterial ${headers} ${sources})

if(USING_BOOST)
	target_link_libraries(TestMaterial pz boost_unit_test_framework)
else()
	target_link_libraries(TestMaterial pz)
endif()

endif()

