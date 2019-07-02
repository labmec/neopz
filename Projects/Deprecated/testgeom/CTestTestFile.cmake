ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (testgeom testgeom_TutorialComp)

#Checking generated files
#SET (filename "malhateste.txt")
#SET (expectedMD5 "9877a0ba3ab80e103f482c6524b173f7")
#add_test (testgeom_malhateste.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
