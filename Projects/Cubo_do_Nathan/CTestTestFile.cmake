ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (Cubo_do_Nathan Cubo_do_Nathan)

#Checking generated files
#SET (filename "Cube.vtk")
#SET (expectedMD5 "64aa9b689c442c5bf3b0a2b428fc0985")
#add_test (Cubo_do_Nathan_Cube.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "malhaPZ1BC.txt")
#SET (expectedMD5 "fdc16e681800611f9730a5ffcf9d231b")
#add_test (Cubo_do_Nathan_malhaPZ1BC.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "cuboPostProc.scal_vec.0.vtk")
#SET (expectedMD5 "35102cc9f2356ec35b643122ff1a45af")
#add_test (Cubo_do_Nathan_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "cuboPostProc.scal_vec.1.vtk")
#SET (expectedMD5 "d799954960d5de1bcf5cf8ad97fec77c")
#add_test (Cubo_do_Nathan_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "cuboPostProc.scal_vec.2.vtk")
#SET (expectedMD5 "e3e4e20057f756973a4db63ed4a9d110")
#add_test (Cubo_do_Nathan_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "cuboPostProc.scal_vec.3.vtk")
#SET (expectedMD5 "26918708ae4a7df72e820fea7f27b109")
#add_test (Cubo_do_Nathan_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "cuboPostProc.scal_vec.4.vtk")
#SET (expectedMD5 "9a4f94f9aa8dd3d0022c735e05059a79")
#add_test (Cubo_do_Nathan_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "cuboPostProc.scal_vec.5.vtk")
#SET (expectedMD5 "34d1fbfb162223e1c6b9d7a537660e7e")
#add_test (Cubo_do_Nathan_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")


#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
