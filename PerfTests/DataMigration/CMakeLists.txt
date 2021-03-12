set(targetName migration-perf)
add_pz_internal_target(
  NAME ${targetName}
  SOURCES main.cpp pzparallel.h)
#linking additional libraries
target_link_libraries(${targetName} PRIVATE Threads::Threads)
