# @file neopz/Projects/Perf/scripts/CMakeList.txt

file(GLOB scripts *.py)

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/test.py"
                "${CMAKE_CURRENT_BINARY_DIR}/test.py")  

set(PERFTEST_SCRIPTS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
  
# Add perf scripts here
add_subdirectory(substruct_tst1)
add_subdirectory(substruct_tst2)
add_subdirectory(substruct_tst3)
add_subdirectory(substruct_tst4)
add_subdirectory(substruct_tst5)