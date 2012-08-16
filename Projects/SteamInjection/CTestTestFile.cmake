ENABLE_TESTING()

#Running the executable and test if it will run ok
ADD_TEST (SteamInjection SteamInjection)

#Checking generated files
#SET (filename "steamInjectionMATHEMATICA.nb")
#SET (expectedMD5 "fe316d1f7271810e0f3806bd2e642aa1")
#add_test (SteamInjection_${filename} "${CMAKE_COMMAND}" "-Dfilename:STRING=${filename}" "-DexpectedMD5:STRING=${expectedMD5}" "-P" "${CMAKE_SOURCE_DIR}/Tests/MD5FileTest.cmake")
#Adding dependency
#set_tests_properties (VisualMatrix_visualmatrix.dx PROPERTIES DEPENDS "VisualMatrix" REQUIRED_FILES "${filename}")
