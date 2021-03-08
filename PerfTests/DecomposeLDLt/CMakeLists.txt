set(targetName Perf-DecomposeLDLt)

add_pz_target(
  NAME ${targetName}
  SOURCES decompose.cpp)
#linking additional libraries
if(USING_LOG4CXX)
  target_link_libraries(${targetName} PRIVATE ${Log4cxx_LIBRARIES})
endif()
