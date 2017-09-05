# @file neopz/UnitTest_PZ/TestMaterial/CMakeLists.txt  -- CMake file for unit test of the matrix module

if (NOT BUILD_COMPLEX_PROJECTS)

    file(GLOB headers *.h)
    file(GLOB sources *.cpp)

    file(COPY "CubeStiffMatrix.txt" DESTINATION ../../../NeoPZ_CMake/UnitTest_PZ/TestMaterial/Debug)
    file(COPY "CubeStiffMatrix.txt" DESTINATION Debug)
    configure_file("CubeStiffMatrix.txt" "CubeStiffMatrix.txt" COPYONLY)

    if(BUILD_UNITTESTING)

            include (CTestTestFile.cmake)

            set(pztestmaterial ${headers} ${sources} PARENT_SCOPE )

            add_executable(TestMaterial ${headers} ${sources})

    if(USING_BOOST)
            IF (WIN32)
                    target_link_libraries(TestMaterial pz ${Boost_LIBRARIES})
            ELSE()
                    target_link_libraries(TestMaterial pz)
            ENDIF()
    else()
            target_link_libraries(TestMaterial pz)
    endif()

    endif()
endif()

