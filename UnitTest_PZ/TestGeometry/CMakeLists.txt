# @file neopz/UnitTest_PZ/TestGeometry/CMakeLists.txt -- CMake file for unit test of geometry module

add_unit_test(TestGeometry
              GeometryUnitTest.cpp
              TPZCurve.cpp
              TPZCurve.h
              TPZSurface.cpp
              TPZSurface.h
              )
