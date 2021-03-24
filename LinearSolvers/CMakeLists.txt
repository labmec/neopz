# @file neopz/LinearSolvers/CMakeLists.txt  -- CMake file for the LinearSolvers module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/LinearSolvers>
                           )

set(headers
    bicg.h
    bicgstab.h
    cg.h
    cgbet.h
    cgs.h
    cheby.h
    gmres.h
    ir.h
    qmr.h
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/LinearSolvers)

target_sources(pz PRIVATE ${headers})
