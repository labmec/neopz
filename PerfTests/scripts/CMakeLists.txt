# @file neopz/Projects/Perf/scripts/CMakeList.txt

file(GLOB scripts *.py)

set(PERFTEST_SCRIPTS_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/test.py"
                "${CMAKE_CURRENT_BINARY_DIR}/test.py")  
  
# Add substruct-perf scripts here
add_subdirectory(substruct_tst01)
add_subdirectory(substruct_tst02)
add_subdirectory(substruct_tst03)
add_subdirectory(substruct_tst04)
add_subdirectory(substruct_tst05)
add_subdirectory(substruct_tst06)
add_subdirectory(substruct_tst07)
add_subdirectory(substruct_tst08)
add_subdirectory(substruct_tst09)
add_subdirectory(substruct_tst10)
add_subdirectory(substruct_tst11)
add_subdirectory(substruct_tst12)
add_subdirectory(substruct_tst13)
add_subdirectory(substruct_tst14)
add_subdirectory(substruct_tst15)
add_subdirectory(substruct_tst16)

# Add skyline-perf scripts here
add_subdirectory(skyline_tst01)
add_subdirectory(skyline_tst02)
add_subdirectory(skyline_tst03)
add_subdirectory(skyline_tst04)
add_subdirectory(skyline_tst05)
add_subdirectory(skyline_tst06)
add_subdirectory(skyline_tst07)
add_subdirectory(skyline_tst08)
add_subdirectory(skyline_tst09)
add_subdirectory(skyline_tst10)
add_subdirectory(skyline_tst11)
add_subdirectory(skyline_tst12)
