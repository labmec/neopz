# @file neopz/Projects/RationalFunction/CMakeLists.txt  -- Simple project to test implementation of rational functions on pyramidal finite element with performance evaluation
set(targetName RationalFunc-Perf)

add_pz_internal_target(
  NAME ${targetName}
  SOURCES main.cpp CedricTest.cpp CedricTest.h)

if(USING_PZ_LOG)
  target_link_libraries(${targetName} PRIVATE ${Log4cxx_LIBRARY})
endif()
