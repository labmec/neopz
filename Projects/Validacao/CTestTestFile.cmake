ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (Validacao Validacao)

#Checking generated files
SET (filename "MalhaDomioTodo.vtk")
SET (expectedMD5 "3e7ff74e99b4066c7be84f7a80fe5f56")
add_test (Validacao_MalhaDomioTodo.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

SET (filename "ErroHdivTodo.txt")
SET (expectedMD5 "75fe6265ffc73c2f0ae5bb94a9793be7")
add_test (Validacao_ErroHdivTodo.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
