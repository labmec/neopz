# @file neopz/UnitTest_PZ/TestMesh/CMakeLists.txt  -- CMake file for unit test of the mesh module
add_unit_test(TestReduced TestReduced.cpp)
# add_executable(TestReduced TestReduced.cpp)
# target_compile_definitions(TestReduced PRIVATE CATCH_CONFIG_MAIN)
# target_link_libraries(TestReduced PUBLIC pz Catch2::Catch2)