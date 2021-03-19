set(targetName Perf-SubStruct)

add_pz_internal_target(
  NAME ${targetName}
  SOURCES substruct.cpp
  FILES cube1.txt 8andares02.txt)

if(USING_PZ_LOG)
  target_link_libraries(${targetName} PRIVATE ${Log4cxx_LIBRARY})
endif()
