ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (CursoPZ1_NumIntegration CursoPZ1_NumIntegration)

#Checking generated files
#SET (filename "PyramidQuad.txt")
#SET (expectedMD5 "8e415a3531100fcd1c1eb7b3c65be5f2")
#add_test (CursoPZClass1_PyramidQuad "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#Adding dependency
#set_tests_properties (CursoPZ_Class1_PyramidQuad.txt PROPERTIES DEPENDS "CursoPZ_Class1" REQUIRED_FILES "${filename}")
