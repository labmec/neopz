ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (QuadraticElements QuadraticElements)

#Checking generated files
#SET (filename "meshAntes.txt")
#SET (expectedMD5 "6f0520a21c89593e50bae04381db91ff")
#add_test (QuadraticElements_meshAntes.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "meshDepois.txt")
#SET (expectedMD5 "d80c65d15bfb6bb1daa899851b034975")
#add_test (QuadraticElements_meshDepois.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")




#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
