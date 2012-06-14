ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (DoubleProjection1D DoubleProjection1D)

#Checking generated files
#SET (filename "cmesh1.txt")
#SET (expectedMD5 "e25f7454012a9b051e9175d1e13ce667")
#add_test (DoubleProjection1D_cmesh1.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "cmesh2.txt")
#SET (expectedMD5 "e25f7454012a9b051e9175d1e13ce667")
#add_test (DoubleProjection1D_cmesh2.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "cmesh12.txt")
#SET (expectedMD5 "f4a2acf1cbc72b8d09c7610c4612d5c7")
#add_test (DoubleProjection1D_cmesh12.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "cmesh22.txt")
#SET (expectedMD5 "af42c902e85b36b659a9746f2b12ac58")
#add_test (DoubleProjection1D_cmesh22.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "gmesh1.txt")
#SET (expectedMD5 "59cf67dce734984841ca22d36d54b0a2")
#add_test (DoubleProjection1D_gmesh1.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "gmesh2.txt")
#SET (expectedMD5 "57ba34cddaab4bd867a9f803ac4f34e7")
#add_test (DoubleProjection1D_gmesh2.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "gmesh3.txt")
#SET (expectedMD5 "2760c687b67c5321b28945059260761c")
#add_test (DoubleProjection1D_gmesh3.txt "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "malhageo1.vtk")
#SET (expectedMD5 "7d244812053420b7ac29c64afc84d7e6")
#add_test (DoubleProjection1D_malhageo1.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "malhageo2.vtk")
#SET (expectedMD5 "7423ba5a9519ff8cc5671297214122a5")
#add_test (DoubleProjection1D_malhageo2.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#SET (filename "malhageoInicial.vtk")
#SET (expectedMD5 "a319762dd4645bdac8d3a1a1d304d37a")
#add_test (DoubleProjection1D_malhageoInicial.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
