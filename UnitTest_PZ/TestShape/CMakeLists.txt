# @file neopz/UnitTest_PZ/TestShape/CMakeLists.txt
# UnitTest of shape functions generation

file(GLOB headers *.h)
file(GLOB sources *.cpp)

if (BUILD_UNITTESTING)

    include(CTestTestFile.cmake)

    set(pztestshape ${headers} ${sources} PARENT_SCOPE)

    add_executable(TestShape ${headers} ${sources})

    
    IF (WIN32)
      target_link_libraries(TestShape pz ${Boost_LIBRARIES})
    ELSE ()
      target_link_libraries(TestShape pz)
      target_link_libraries(TestShape ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    ENDIF ()

endif ()
