# @file neopz/UnitTest_PZ/TestMesh/CMakeLists.txt  -- CMake file for unit test of the mesh module

file(GLOB headers *.h)
file(GLOB sources1 TestMHM.cpp TPZMatLaplacianHybrid.cpp)


if(BUILD_UNITTESTING)

	include (CTestTestFile.cmake)

	set(pztestmmhm ${headers} ${sources1} PARENT_SCOPE)

        if(NOT BUILD_COMPLEX_PROJECTS)
	    add_executable(TestMMHM ${headers} ${sources1})
        endif()

if(USING_BOOST)
	IF (WIN32)
            if(NOT BUILD_COMPLEX_PROJECTS)
		target_link_libraries(TestMMHM pz ${Boost_LIBRARIES})
            endif()
	ELSE()
            if(NOT BUILD_COMPLEX_PROJECTS)
		target_link_libraries(TestMMHM pz)
            endif()
	ENDIF()
else()
    if(NOT BUILD_COMPLEX_PROJECTS)
	target_link_libraries(TestMHM pz)
    endif()
endif()

endif()

