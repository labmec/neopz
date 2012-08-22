ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (testshape_Tutorial testshape_Tutorial)

#Checking generated files
#SET (filename "malha.txt")
#SET (expectedMD5 "255013f097df2aa8faac7d42f204043d")
#add_test (testshape_Tutorial_malha.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
