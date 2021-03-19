set(targetName Perf-NumaTst)

add_pz_internal_target(
  NAME ${targetName}
  SOURCES numatst.cpp)
#linking additional libraries
target_link_libraries(${targetName} PRIVATE Threads::Threads)
if(USING_PZ_LOG)
  target_link_libraries(${targetName} PRIVATE ${Log4cxx_LIBRARY})
endif()
