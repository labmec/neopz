# @file neopz/UnitTest_PZ/TestMesh/CMakeLists.txt  -- CMake file for unit test of the mesh module

add_unit_test(TestHDivConstant TestHDivConstant.cpp TPZMatL2Product.cpp)
add_unit_test(TestHDiv TestHDiv.cpp)
# add_unit_test(TestCondensedEl TestCondensedElement.cpp) NEEDS FIX
add_unit_test(TestErrorAnalysis TestError.cpp)
add_unit_test(TestKernelHDiv TestKernelHDiv.cpp)
add_unit_test(TestHangingNode TestHangingNode.cpp)
add_unit_test(TestSideShapeF TestSideShapeF.cpp)