# @file neopz/UnitTest_PZ/TestMatrix/CMakeLists.txt  -- CMake file for unit test of tensor
add_unit_test(TestTensor TestTensor.cpp)
configure_file(tensor_and_eigensystem.txt tensor_and_eigensystem.txt COPYONLY)