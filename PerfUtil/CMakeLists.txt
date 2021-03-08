# @file neopz/PerfUtil/CMakeLists.txt  -- CMake file for the PerfUtil module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/PerfUtil>
                          )

set(headers
        arglib.h
        clock_timer.h
        csvtable.h
        cycle_timer.h
        pz_gettime.h
        run_stats_table.h
        stats_recorder.h
        timing_analysis.h
        )

set(sources arglib.cpp)

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/PerfUtil)
target_sources(pz PRIVATE ${headers} ${sources})
