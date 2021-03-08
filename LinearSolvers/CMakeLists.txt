# @file neopz/LinearSolvers/CMakeLists.txt  -- CMake file for the LinearSolvers module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/LinearSolvers>
                          )

set(headers
    TPZPardisoHandler.h
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

set(sources TPZPardisoHandler.cpp)

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/LinearSolvers)

target_sources(pz PRIVATE ${headers} ${sources})
