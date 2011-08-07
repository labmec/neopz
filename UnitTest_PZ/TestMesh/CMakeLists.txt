# @file neopz/UnitTest_PZ/TestMesh/CMakeLists.txt  -- CMake file for unit test of the mesh module

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

	set(pztestmesh ${headers} ${sources} PARENT_SCOPE)

	add_executable(TestMesh ${headers} ${sources})

	target_link_libraries(TestMesh pz)

endif()

