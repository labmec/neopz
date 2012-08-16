ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (MultiphysicsSimulation MultiphysicsSimulation)

#Checking generated files
#SET (filename "Solution.out")
#SET (expectedMD5 "abd7ef8ed016c47d0ed26d76ebda2544")
#add_test (MultiphysicsSimulation_Solution.out "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
