ENABLE_TESTING()

#Running the executable and test if it will run ok
#ADD_TEST (Plasticity Plasticity)

#Checking generated files
#SET (filename "barmesh.txt")
#SET (expectedMD5 "0d5e7dd23959b4ddcff549ed9fff9ec7")
#add_test (Plasticity_barmesh.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "barmesh.vtk")
#SET (expectedMD5 "03a6c0afb084ca317134b16950e3473b")
#add_test (Plasticity_barmesh.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "hyper.scal_vec.0.vtk")
#SET (expectedMD5 "7e810b6f8e57a812d449caf30ffc101b")
#add_test (Plasticity_hyper.scal_vec.0.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "hyper.txt")
#SET (expectedMD5 "c15108fe29e5beadc2fea53882b5638b")
#add_test (Plasticity_hyper.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")


#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
