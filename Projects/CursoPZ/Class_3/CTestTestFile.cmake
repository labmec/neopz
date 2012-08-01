ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (CursoPZ3_Jacobian CursoPZ3_Jacobian)

#Checking generated files
#SET (filename "visualmatrix.dx")
#SET (expectedMD5 "654b8ae8045e237d00a60ccbc7032b1c")
#add_test (Vi89459489495x "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
