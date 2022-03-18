# @file neopz/UnitTest_PZ/TestMesh/CMakeLists.txt  -- CMake file for unit test of the mesh module

add_unit_test(TestHDivConstant TestHDivConstant.cpp)
add_unit_test(TestHDiv TestHDiv.cpp)
# add_unit_test(TestCondensedEl TestCondensedElement.cpp)
add_unit_test(TestErrorAnalysis TestError.cpp)
add_unit_test(TestKernelHDiv TestKernelHDiv.cpp)