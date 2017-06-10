# @file neopz/UnitTest_PZ/TestMesh/CMakeLists.txt  -- CMake file for unit test of the mesh module

file(GLOB headers *.h)
file(GLOB sources1 TestHDiv.cpp)
file(GLOB sources2 TestCondensedElement.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestmesh1 ${headers} ${sources1} PARENT_SCOPE)
	set(pztestmesh2 ${headers} ${sources2} PARENT_SCOPE)

        if(NOT BUILD_COMPLEX_PROJECTS)
	    add_executable(TestMesh1 ${headers} ${sources1})
        endif()
	add_executable(TestMesh2 ${headers} ${sources2})

if(USING_BOOST)
	IF (WIN32)
            if(NOT BUILD_COMPLEX_PROJECTS)
		target_link_libraries(TestMesh1 pz ${Boost_LIBRARIES})
            endif()
		target_link_libraries(TestMesh2 pz ${Boost_LIBRARIES})
	ELSE()
            if(NOT BUILD_COMPLEX_PROJECTS)
		target_link_libraries(TestMesh1 pz)
            endif()
		target_link_libraries(TestMesh2 pz)
	ENDIF()
else()
    if(NOT BUILD_COMPLEX_PROJECTS)
	target_link_libraries(TestMesh1 pz)
    endif()
    target_link_libraries(TestMesh2 pz)
endif()

endif()

