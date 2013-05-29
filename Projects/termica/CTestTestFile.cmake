ENABLE_TESTING()

#Running the executable and test if it will run ok
#ADD_TEST (termica termica)

#Checking generated files
#SET (filename matrixstruct.vtk)
#SET (expectedMD5 edd29a9dc8afee12f066590186c3cb9d)
#add_test (termica_matrixstruct.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#Adding dependency
#set_tests_properties (termica_matrixstruct.vtk PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
#ADD_TEST(termica_matrixstruct.vtk "/Applications/CMake 2.8-5.app/Contents/bin/cmake" "-Dfilename:STRING=matrixstruct.vtk" "-DexpectedMD5:STRING=edd29a9dc8afee12f066590186c3cb9d" "-P" "/Users/emerson/Documents/pz/Tests/MD5FileTest.cmake")
