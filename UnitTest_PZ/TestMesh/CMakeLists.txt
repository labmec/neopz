# @file neopz/UnitTest_PZ/TestMesh/CMakeLists.txt  -- CMake file for unit test of the mesh module

if (NOT BUILD_COMPLEX_PROJECTS)

    set(TestSources1 TestHDiv.cpp)
    set(TestSources2 TestCondensedElement.cpp)

    # add_unit_test(TestMesh1 ${TestSources1})
    # add_unit_test(TestMesh2 ${TestSources2})

endif ()
