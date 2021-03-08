set(targetName Perf-SubStruct)

add_pz_target(
  NAME ${targetName}
  SOURCES substruct.cpp
  FILES cube1.txt 8andares02.txt)

if(USING_LOG4CXX)
  target_link_libraries(${targetName} PRIVATE ${Log4cxx_LIBRARY})
endif()
