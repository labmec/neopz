# @file neopz/Solvers/EigenSolvers/CMakeLists.txt  -- CMake file for the EigenSolvers module

target_include_directories(pz PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Solvers/EigenSolvers>
    )

set(public_headers
    TPZLapackEigenSolver.h
    TPZKrylovEigenSolver.h
    TPZSpectralTransform.h
    )
set(headers
    ${public_headers}
    )

set(sources
    TPZLapackEigenSolver.cpp
    TPZKrylovEigenSolver.cpp
    TPZSpectralTransform.cpp
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Solvers/EigenSolvers)

target_sources(pz PRIVATE ${headers} ${sources})
#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)