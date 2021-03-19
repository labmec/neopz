########### next target ###############
set(targetName ParallelPerformance)
add_pz_internal_target(
  NAME ${targetName}
  SOURCES main.cpp)
#linking additional libraries
if(USING_PZ_LOG)
  target_link_libraries(${targetName} PRIVATE ${Log4cxx_LIBRARY})
 endif()
