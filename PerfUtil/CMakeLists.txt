# @file neopz/PerfUtil/CMakeLists.txt  -- CMake file for the PerfUtil module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           #$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/PerfUtil>
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

  set(sources
    arglib.cpp
    stats_recorder.cpp
    )


#@orlandini: until we make sure that they are actually still useful i guess we dont need to install them  
# install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PerfUtil)

target_sources(pz PRIVATE ${headers} ${sources})
