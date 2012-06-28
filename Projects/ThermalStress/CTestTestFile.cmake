ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (ThermalStress ThermalStress)

#Checking generated files
#SET (filename "perkins.scal_vec.0.vtk")
#SET (expectedMD5 "6c2523fef65390fdf2196afe09d09bff")
#add_test (ThermalStress_perkins.scal_vec.0.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
