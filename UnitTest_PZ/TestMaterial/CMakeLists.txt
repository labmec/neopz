# @file neopz/UnitTest_PZ/TestMaterial/CMakeLists.txt  -- CMake file for unit test of the material module

if (NOT BUILD_COMPLEX_PROJECTS)
    configure_file("CubeStiffMatrix.txt" "CubeStiffMatrix.txt" COPYONLY)
    add_unit_test(TestMaterial TestMaterial.cpp)
endif ()
