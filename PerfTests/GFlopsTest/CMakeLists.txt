set(targetName GFlopsTest)
add_pz_target(
  NAME ${targetName}
  SOURCES gflopstest.cpp)
#linking additional libraries
target_link_libraries(${targetName} PRIVATE Threads::Threads)

