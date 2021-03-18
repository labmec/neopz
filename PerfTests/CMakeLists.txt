# @file neopz/PerfTests/CMakeList.txt

if(BUILD_PERF_TESTS)

  if(NOT UNIX)  
     # UNIX is TRUE on all UNIX-like OS's, including Apple OS X and CygWin
     message (FATAL_ERROR "ERROR: PERFTESTS may only be used on UNIX-like OS's")
  endif()

  # Search for the directory containing small data sets. Small data sets can be
  # kept together with the source code, hence ./small_data

  find_path(PERFTEST_SMALL_DATA_DIR libpz_small_datadir.readme PATHS ./small_data)
  if(NOT PERFTEST_SMALL_DATA_DIR)
  	message (FATAL_ERROR "ERROR: Please specify a valid directory for PERFTEST_SMALL_DATA_DIR")
  endif()

  # Search for the directory containing large data sets
  find_path(PERFTEST_LARGE_DATA_DIR libpz_large_datadir.readme PATHS ${EXTRA_SEARCH_DIRS})
  find_path(PERFTEST_LARGE_DATA_DIR libpz_large_datadir.readme PATHS ../../neopz_perftest_large_data)
  find_path(PERFTEST_LARGE_DATA_DIR libpz_large_datadir.readme PATHS ../neopz_perftest_large_data)

  if(NOT PERFTEST_LARGE_DATA_DIR)
  	#message (WARNING "WARNING: Could not find the path for large input and reference data sets "
	#  "PERFTEST_LARGE_DATA_DIR. Performance Tests that rely on large data sets may not work!")
  endif()

  # This variable is used inside the run scripts
  set(PERFTEST_BASE_DIR ${CMAKE_CURRENT_BINARY_DIR})

  # Include the subdirectories
  # add_subdirectory(progs)
  add_subdirectory(scripts)
  
  message("PERFTEST_BASE_DIR = ${PERFTEST_APPS_DIR}")
  message("PERFTEST_SMALL_DATA_DIR = ${PERFTEST_SMALL_DATA_DIR}")
  message("PERFTEST_LARGE_DATA_DIR = ${PERFTEST_LARGE_DATA_DIR}")

  # Old structure -- TODO: Move these to progs
  # Add perf test applications here
  if(USING_METIS)
    add_subdirectory(SubStruct)
  endif()
  add_subdirectory(DecomposeLDLt)
  add_subdirectory(NumaTests)
  add_subdirectory(GFlopsTest)
  if(USING_TBB)
    add_subdirectory(TBBTests)
  endif()
  add_subdirectory(RationalFunc)
  add_subdirectory(DataMigration)
  add_subdirectory(ParallelAssembly)

  # Add your run scripts into RunScripts/CMakeLists.txt
  add_subdirectory(RunScripts)

endif(BUILD_PERF_TESTS)
		
