# @file neopz/UnitTest_PZ/TestMultiThreading/CMakeLists.txt  -- CMake file for unit test of multithreading in NeoPZ

file(GLOB headers *.h)
file(GLOB sources *.cpp)


if(BUILD_UNITTESTING)

  add_test(TestMultithreading TestMultithreading)
  
  add_executable(TestMultithreading ${headers} ${sources})
  if (WIN32)
    target_link_libraries(TestMultithreading pz ${Boost_LIBRARIES})
  else()
    target_link_libraries(TestMultithreading pz)
    target_link_libraries(TestMultithreading ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  endif()
endif(BUILD_UNITTESTING)
