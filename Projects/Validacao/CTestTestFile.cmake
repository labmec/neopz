ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (Validacao Validacao)

#Checking generated files
#SET (filename "MalhaDomioTodo.vtk")
#SET (expectedMD5 "3e7ff74e99b4066c7be84f7a80fe5f56")
#add_test (Validacao_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "ErroHdivTodo.txt")
#SET (expectedMD5 "85551d58936589043028c3b440b36df6")
#add_test (Validacao_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
  
#SET (filename "GraficoH1Todo.scal_vec.0.vtk")
#SET (expectedMD5 "743f3c285f2c9baf6d31606bf4d3fb48")
#add_test (Validacao_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
