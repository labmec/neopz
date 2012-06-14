ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (SteamInjection SteamInjection)

#Checking generated files
SET (filename "teamInjectionMATHEMATICA.nb")
SET (expectedMD5 "fe9ee68844ff96272fba5186e8207c41")
add_test (SteamInjection_teamInjectionMATHEMATICA.nb "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
