# @file neopz/UnitTest_PZ/TestMaterial/CMakeLists.txt  -- CMake file for unit test of the material module

configure_file("CubeStiffMatrix.txt" "CubeStiffMatrix.txt" COPYONLY)
add_unit_test(TestMaterial TestMaterial.cpp)

add_unit_test(TestDarcyFlow TestDarcyFlow.cpp)