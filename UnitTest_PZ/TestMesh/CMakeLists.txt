# @file neopz/UnitTest_PZ/TestMesh/CMakeLists.txt  -- CMake file for unit test of the mesh module

if (NOT BUILD_COMPLEX_PROJECTS)
    # add_unit_test(TestHDiv TestHDiv.cpp)
    #add_unit_test(TestCondensedEl TestCondensedElement.cpp)
    add_unit_test(TestErrorAnalysis TestError.cpp)
endif ()
