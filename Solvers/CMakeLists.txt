# @file neopz/Solvers/CMakeLists.txt  -- CMake file for Solvers module. 

target_include_directories(pz PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Solvers>
    )
set(public_headers
    TPZSolver.h
    TPZMatrixSolver.h
    TPZEigenSort.h
    TPZEigenSolver.h
    TPZPardisoSolver.h
    TPZIdentitySolver.h
    TPZJacobiPrecond.h
    pzstepsolver.h
    pzseqsolver.h
    )

set(sources
    TPZSolver.cpp
    TPZMatrixSolver.cpp
    TPZEigenSolver.cpp
    TPZPardisoSolver.cpp
    TPZJacobiPrecond.cpp
    pzstepsolver.cpp
    pzseqsolver.cpp
    )

if(USING_MKL)
    list(APPEND headers
      TPZPardisoSolver.h
    )
    list(APPEND sources
      TPZPardisoSolver.cpp
      )
endif()
    
install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Solvers)

target_sources(pz PRIVATE ${public_headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

add_subdirectory(LinearSolvers)
add_subdirectory(EigenSolvers)
add_subdirectory(Multigrid)
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)