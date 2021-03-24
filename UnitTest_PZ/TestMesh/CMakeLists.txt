# @file neopz/UnitTest_PZ/TestMesh/CMakeLists.txt  -- CMake file for unit test of the mesh module

if (NOT BUILD_COMPLEX_PROJECTS)
    add_unit_test(TestMesh1 TestHDiv.cpp)
    add_unit_test(TestMesh2 TestCondensedElement.cpp)
endif ()
