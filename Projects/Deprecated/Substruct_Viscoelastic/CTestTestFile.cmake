ENABLE_TESTING()

#Running the executable and test if it will run ok
#ADD_TEST (SubStruct_Viscoelastic Substruct_Viscoelastic)

#Checking generated files
#SET (filename "CoarseMatrix.vtk")
#SET (expectedMD5 "21b4804f7ffac977c2ba8c25e699f3c7")
#add_test (SubStruct_Viscoelastic_CoarseMatrix.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "Cube.vtk")
#SET (expectedMD5 "64aa9b689c442c5bf3b0a2b428fc0985")
#add_test (SubStruct_Viscoelastic_Cube.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "dohrmann_visco.scal_vec.0.vtk")
#SET (expectedMD5 "c1eea00cf5dbcf7848104ef11a9c187a")
#add_test (SubStruct_Viscoelastic_dohrmann_visco.scal_vec.0.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "malhaPZ1BC.txt")
#SET (expectedMD5 "fdc16e681800611f9730a5ffcf9d231b")
#add_test (SubStruct_Viscoelastic_malhaPZ1BC.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "partition.vtk")
#SET (expectedMD5 "81d786560029ef762e138d658157da10")
#add_test (SubStruct_Viscoelastic_partition.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "partitionbefore.vtk")
#SET (expectedMD5 "85012abb71d8463225273eeebe01ee7e")
#add_test (SubStruct_Viscoelastic_partitionbefore.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "SubMatrix82.vtk")
#SET (expectedMD5 "12257ed7f55651db4d0a3ab8abd8da80")
#add_test (SubStruct_Viscoelastic_SubMatrix82.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "SubMatrix83.vtk")
#SET (expectedMD5 "756d16950e44e37eff9406c0aa7a0549")
#add_test (SubStruct_Viscoelastic_SubMatrix83.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "SubMatrix84.vtk")
#SET (expectedMD5 "bf40abb00e996049ddd33e7df46f71bc")
#add_test (SubStruct_SubMatrix84.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "SubMatrixInternal82.vtk")
#SET (expectedMD5 "8b15752850a70992f6a41a525c4d5fe7")
#add_test (SubStruct_Viscoelastic_SubMatrixInternal82.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "SubMatrixInternal83.vtk")
#SET (expectedMD5 "cc74a68792733b3b1ffdf692c5068380")
#add_test (SubStruct_Viscoelastic_SubMatrixInternal83.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "SubMatrixInternal84.vtk")
#SET (expectedMD5 "5e755002bf76111e624f1c1677d27d41")
#add_test (SubStruct_Viscoelastic_SubMatrixInternal84.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "dohrmann_visco.scal_vec.1.vtk")
#SET (expectedMD5 "d5e4dba3d8759e2d76bb3b7900f10478")
#add_test (SubStruct_Viscoelastic_dohrmann_visco.scal_vec.1.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "dohrmann_visco.scal_vec.2.vtk")
#SET (expectedMD5 "a7e05b789eaac4cb67d44703739e8176")
#add_test (SubStruct_Viscoelastic_dohrmann_visco.scal_vec.2.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "dohrmann_visco.scal_vec.3.vtk")
#SET (expectedMD5 "a58e636f239598831518bb8d000acc53")
#add_test (SubStruct_Viscoelastic_dohrmann_visco.scal_vec.3.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "dohrmann_visco.scal_vec.4.vtk")
#SET (expectedMD5 "d77a2aa2836bf580b9c48a9b619bcb13")
#add_test (SubStruct_Viscoelastic_dohrmann_visco.scal_vec.4.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "PointMesh.vtk")
#SET (expectedMD5 "8141dfaa69e81acd540f63d959d9e230")
#add_test (SubStruct_Viscoelastic_PointMesh.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")


#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
