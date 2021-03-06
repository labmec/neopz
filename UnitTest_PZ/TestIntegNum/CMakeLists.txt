# @file neopz/UnitTest_PZ/TestIntegNum/CMakeLists.txt  -- CMake file for unit test of numerical integration

set(TestSources
        IntegNumUnitTest.cpp
        GeometryUnitTest.cpp
        TPZCurve.cpp
        TPZCurve.h
        )
add_unit_test(TestIntegNum ${TestSources})
