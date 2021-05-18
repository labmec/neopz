# @file neopz/Solvers/CMakeLists.txt  -- CMake file for Solvers module. 

target_include_directories(pz PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Solvers>
    )
set(public_headers
    TPZSolver.h
    TPZMatrixSolver.h
    TPZEigenSolver.h
    TPZPardisoSolver.h
    pzstepsolver.h
    pzseqsolver.h
    )
set(headers
    ${public_headers}
    TPZCopySolve.h
    )

set(sources
    TPZSolver.cpp
    TPZMatrixSolver.cpp
    TPZEigenSolver.cpp
    TPZPardisoSolver.cpp
    pzstepsolver.cpp
    pzseqsolver.cpp
    )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Solvers)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

add_subdirectory(LinearSolvers)
add_subdirectory(EigenSolvers)
add_subdirectory(Multigrid)
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)