# @file neopz/UnitTest_PZ/TestGeometry/CMakeLists.txt -- CMake file for unit test of geometry module

if (USING_FAD)
    add_unit_test(TestGeometry
            GeometryUnitTest.cpp
            TPZCurve.cpp
            TPZCurve.h
            TPZSurface.cpp
            TPZSurface.h
            )
endif ()
