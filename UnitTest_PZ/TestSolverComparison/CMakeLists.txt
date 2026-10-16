# @file neopz/UnitTest_PZ/TestSolverComparison/CMakeLists.txt

if(USING_MKL AND USING_MUMPS)
    add_unit_test(TestSolverComparison TestSolverComparison.cpp)
endif()
