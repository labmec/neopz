# @file neopz/UnitTest_PZ/TestMesh/CMakeLists.txt  -- CMake file for unit test of the mesh module

file(GLOB headers *.h)
file(GLOB sources1 TestHDiv.cpp)
file(GLOB sources2 TestCondensedElement.cpp)


if(BUILD_UNITTESTING AND NOT BUILD_COMPLEX_PROJECTS)

	include (CTestTestFile.cmake)

	set(pztestmesh1 ${headers} ${sources1} PARENT_SCOPE)
	set(pztestmesh2 ${headers} ${sources2} PARENT_SCOPE)

        if(NOT BUILD_COMPLEX_PROJECTS)
	    add_executable(TestMesh1 ${headers} ${sources1})
        endif()
	add_executable(TestMesh2 ${headers} ${sources2})
        if (WIN32)
	  target_link_libraries(TestMesh1 pz ${Boost_LIBRARIES})
          target_link_libraries(TestMesh2 pz ${Boost_LIBRARIES})
	else()
	  target_link_libraries(TestMesh1 pz)
          target_link_libraries(TestMesh1 ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
          target_link_libraries(TestMesh2 pz)
          target_link_libraries(TestMesh2 ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
        endif()

endif()

