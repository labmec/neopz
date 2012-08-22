ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (Adapt_HP_Jorge Adapt_HP_Jorge)

#Checking generated files
#SET (filename "aftercmesh.vtk")
#SET (expectedMD5 "fe1fc811cb3e6af9b4aeb794796ab62c")
#add_test (Adapt_HP_Jorge_aftercmesh.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "malhas.vtk")
#SET (expectedMD5 "ac852a8c5fb64cf131172ce76f019a05")
#add_test (Adapt_HP_Jorge_malhas.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "meshes.vtk")
#SET (expectedMD5 "6deb62bb33c0fd2a7cf2ea35eed64a69")
#add_test (Adapt_HP_Jorge_meshes.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "meshextruded.vtk")
#SET (expectedMD5 "df0c3b5c837bb96115dae7433037788c")
#add_test (Adapt_HP_Jorge_meshextruded.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "meshInitialSol.vtk")
#SET (expectedMD5 "c79fb2cf3262c0040e9f61d84aec2c96")
#add_test (Adapt_HP_Jorge_meshInitialSol.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "meshrefined.vtk")
#SET (expectedMD5 "27cac4fb3f9bce8fbed7da6731159d1c")
#add_test (Adapt_HP_Jorge_meshrefined.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "original.vtk")
#SET (expectedMD5 "28ec16299923d1708a80356d655d4899")
#add_test (Adapt_HP_Jorge_original.vtk "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#SET (filename "testeinicialadapt_0.dx")
#SET (expectedMD5 "c2a5df37dfcc933c8e8b576f0e6f444e")
#add_test (Adapt_HP_Jorge_testeinicialadapt_0.dx "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")

#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
