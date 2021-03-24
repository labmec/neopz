# @file neopz/UnitTest_PZ/CMakeList.txt  -- First Cmake file for Unit Testing modules

enable_testing()

add_subdirectory(TestFad)
add_subdirectory(TestHCurl)
add_subdirectory(TestTopology)
add_subdirectory(TestMatrix)
#add_subdirectory(TestMesh)
add_subdirectory(TestMaterial)
add_subdirectory(TestIntegNum)
add_subdirectory(TestIterator)
add_subdirectory(TestGeometry)
add_subdirectory(TestPlasticity)
add_subdirectory(TestTensor)
add_subdirectory(TestRandom)
add_subdirectory(TestBlend)
add_subdirectory(TestShape)
add_subdirectory(TestMultithreading)
