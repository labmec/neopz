ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (CursoPZ_Class5 CursoPZ_Class5)

#Checking generated files
SET (filename "all.dat")
SET (expectedMD5 "d41d8cd98f00b204e9800998ecf8427e")
add_test (CursoPZ_Class5_all.dat "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
