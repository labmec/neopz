# @file neopz/Analysis/CMakeLists.txt  -- CMake file for Analysis module. 

target_include_directories(pz PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Analysis>
    )
set(public_headers
    TPZAnalysis.h
    TPZLinearAnalysis.h
    TPZEigenAnalysis.h
    pznonlinanalysis.h
    pzsmanal.h #TPZSubMeshAnalysis
  )

set(headers
    ${public_headers}
    TPZNLMultGridAnalysis.h
    pzmganalysis.h
    pztransientanalysis.h
    pzeuleranalysis.h
    pzsmfrontalanal.h
    )

set(sources
    TPZAnalysis.cpp
    TPZLinearAnalysis.cpp
    TPZEigenAnalysis.cpp
    pznonlinanalysis.cpp
    pzsmanal.cpp
    TPZNLMultGridAnalysis.cpp
    pzmganalysis.cpp
    pztransientanalysis.cpp
    pzeuleranalysis.cpp
    pzsmfrontalanal.cpp
    )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Analysis)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

#black oil
add_subdirectory(BlackOil)
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)