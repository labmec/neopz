set(targetName Perf-DecomposeLDLt)

add_pz_internal_target(
  NAME ${targetName}
  SOURCES decompose.cpp)
#linking additional libraries
if(USING_PZ_LOG)
  target_link_libraries(${targetName} PRIVATE ${Log4cxx_LIBRARIES})
endif()
