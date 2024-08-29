# @file neopz/UnitTest_PZ/TestMatrix/CMakeLists.txt  -- CMake file for unit test of the matrix module

add_unit_test(TestMatrix TestMatrix.cpp)
add_unit_test(TestMatrix_Self MatrixSelfTest.cpp TestMatrixHeaders.h)
add_unit_test(TestMatrix_Algebra MatrixAlgebraTest.cpp TestMatrixHeaders.h)
if(USING_LAPACK)
    add_unit_test(TestMatrix_Eigen MatrixEigenTest.cpp TestMatrixHeaders.h)
endif()
add_unit_test(TestMatrix_Window MatrixWindowTest.cpp)
