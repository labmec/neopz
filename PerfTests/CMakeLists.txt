# @file neopz/PerfTests/CMakeList.txt
  
if(BUILD_PERF_TESTS)

  if(NOT UNIX)  
     # UNIX is TRUE on all UNIX-like OS's, including Apple OS X and CygWin
     message (FATAL_ERROR "ERROR: PERFTESTS may only be used on UNIX-like OS's")
  endif()

  find_path(PERFTEST_DATA_DIR libpz_datadir.readme PATHS ${SEARCH_DIRS})
  find_path(PERFTEST_DATA_DIR libpz_datadir.readme PATHS ../../neopz_perf_data)
  find_path(PERFTEST_DATA_DIR libpz_datadir.readme PATHS ./neopz_perf_data)

  if(NOT PERFTEST_DATA_DIR)
  	message (FATAL_ERROR "ERROR: Please specify a valid directory for PERFTEST_DATA_DIR")
  endif()

  # This variable is used inside the run scripts
  set(PERFTEST_APPS_DIR ${CMAKE_CURRENT_BINARY_DIR})

  # Add perf test applications here
  add_subdirectory(SubStruct)
  add_subdirectory(GFlopsTest)

  # Add your run scripts into RunScripts/CMakeLists.txt
  add_subdirectory(RunScripts)

  message("PERFTEST_APPS_DIR = ${PERFTEST_APPS_DIR}")
  message("PERFTEST_DATA_DIR = ${PERFTEST_DATA_DIR}")

endif(BUILD_PERF_TESTS)
		
