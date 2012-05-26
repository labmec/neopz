# @file neopz/PerfTests/CMakeList.txt
  
if(BUILD_PERF_TESTS)

  # This variable is used inside the run scripts
  set(PERFTEST_APPS_DIR ${CMAKE_CURRENT_BINARY_DIR})

  # Add perf test applications here
  add_subdirectory(SubStruct)

  # Add your run scripts into RunScripts/CMakeLists.txt
  add_subdirectory(RunScripts)

  message("PERFTEST_APPS_DIR = ${PERFTEST_APPS_DIR}")
  message("PERFTEST_DATA_DIR = ${PERFTEST_DATA_DIR}")

endif(BUILD_PERF_TESTS)
		
