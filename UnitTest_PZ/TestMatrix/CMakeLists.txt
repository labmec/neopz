# @file neopz/UnitTest_PZ/TestMatrix/CMakeLists.txt  -- CMake file for unit test of the matrix module

add_unit_test(TestMatrix TestMatrix.cpp)
add_unit_test(TestMatrix_Self MatrixSelfTest.cpp)
# add_unit_test(TestMatrix_Algebra MatrixAlgebraTest.cpp)
# add_unit_test(TestMatrix_Eigen MatrixEigenTest.cpp)
# add_unit_test(TestMatrix_Windowing MatrixWindowingTest.cpp)
