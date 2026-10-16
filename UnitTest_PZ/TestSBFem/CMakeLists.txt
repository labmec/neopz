# @file neopz/UnitTest_PZ/TestMultiThreading/CMakeLists.txt  -- CMake file for unit test of multithreading in NeoPZ

add_unit_test(TestSBFem SBFemUnitTest.cpp)

add_unit_test(TestSBFemHdiv
    SBFemHdivUnitTest.cpp
    #materials for SBFemHdiv
    TPZDarcySBFemHdiv.cpp
)