# @file neopz/UnitTest_PZ/TestIntegNum/CMakeLists.txt  -- CMake file for unit test of numerical integration

add_unit_test(TestIntegNum
        IntegNumUnitTest.cpp
        TPZCurve.cpp
        TPZCurve.h
        )

configure_file(Line.txt ${CMAKE_CURRENT_BINARY_DIR}/Line.txt COPYONLY)
configure_file(Triangle.txt ${CMAKE_CURRENT_BINARY_DIR}/Triangle.txt COPYONLY)
configure_file(NumInteBeingChecked.nb ${CMAKE_CURRENT_BINARY_DIR}/NumInteBeingChecked COPYONLY)
