# @file neopz/Matrix/CMakeLists.txt  -- CMake file for the Matrix module

target_include_directories(pz PUBLIC 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Matrix>
    )

set(headers
    PZMatrixMarket.h
    doxmatrix.h
    pzbndmat.h
    pzfmatrix.h
    TPZFMatrixRef.h
    pzmatrix.h
    pzbasematrix.h
    pzskylmat.h
    pzspblockdiagpivot.h
    pzsysmp.h
    tpzsparseblockdiagonal.h
    pzblock.h
    pzdiffmatrix.h
    pzsbndmat.h
    pzsfulmat.h
    pzskylnsymmat.h
    pzstencil.h
    pztrnsform.h
    tpzverysparsematrix.h
    pzblockdiag.h
    pzmatred.h
    pzshtmat.h
    pzysmp.h
    TPZSolutionMatrix.h
    TPZTensor.h
    TPZTensorTranslator.h
    )

set(sources
    PZMatrixMarket.cpp
    pzblockdiag.cpp
    pzfmatrix.cpp
    pzmatrix.cpp
    pzbasematrix.cpp
    pzskylmat.cpp
    pzspblockdiagpivot.cpp
    pzsysmp.cpp
    tpzsparseblockdiagonal.cpp
    pzbndmat.cpp
    pzsbndmat.cpp
    pzsfulmat.cpp
    pzskylnsymmat.cpp
    pzstencil.cpp
    pztrnsform.cpp
    tpzverysparsematrix.cpp
    pzblock.cpp
    pzmatred.cpp
    pzshtmat.cpp
    pzysmp.cpp
    TPZSolutionMatrix.cpp
    TPZTensor.cpp
    TPZTensorTranslator.cpp
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Matrix)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)