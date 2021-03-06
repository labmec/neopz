# @file neopz/UnitTest_PZ/TestGeometry/CMakeLists.txt -- CMake file for unit test of geometry module

if (USING_FAD)
    set(TestSources
            GeometryUnitTest.cpp
            TPZCurve.cpp
            TPZCurve.h
            TPZSurface.cpp
            TPZSurface.h
            )
    add_unit_test(TestGeometry ${TestSources})
endif ()
