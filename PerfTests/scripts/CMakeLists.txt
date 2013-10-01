# @file neopz/Projects/Perf/scripts/CMakeList.txt

file(GLOB scripts *.py)

set(PERFTEST_SCRIPTS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/test.py"
                "${CMAKE_CURRENT_BINARY_DIR}/test.py")  
  
# Add substruct-perf scripts here
add_subdirectory(substruct_tst1)
add_subdirectory(substruct_tst2)
add_subdirectory(substruct_tst3)
add_subdirectory(substruct_tst4)
add_subdirectory(substruct_tst5)
add_subdirectory(substruct_tst06)
add_subdirectory(substruct_tst07)
add_subdirectory(substruct_tst08)
add_subdirectory(substruct_tst09)
add_subdirectory(substruct_tst10)
add_subdirectory(substruct_tst11)
add_subdirectory(substruct_tst12)
add_subdirectory(substruct_tst13)
add_subdirectory(substruct_tst14)

# Add skyline-perf scripts here
add_subdirectory(skyline_tst1)
add_subdirectory(skyline_tst2)
add_subdirectory(skyline_tst3)
add_subdirectory(skyline_tst4)
add_subdirectory(skyline_tst5)
add_subdirectory(skyline_tst7)
add_subdirectory(skyline_tst8)
