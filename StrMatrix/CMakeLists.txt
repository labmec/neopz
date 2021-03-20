# @file neopz/StrMatrix/CMakeLists.txt  -- CMake file for the StrMatrix module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/StrMatrix>
                           )
#@TODO: check if these public headers are exaggerated
set(public_headers
    TPZStructMatrixBase.h
    #parallel scheme
    pzstrmatrixot.h
    pzstrmatrix.h
    pzstrmatrixflowtbb.h
    #matrix storage format
    pzbstrmatrix.h
    TPZSpStructMatrix.h
    pzbdstrmatrix.h
    pzsbstrmatrix.h
    pzskylstrmatrix.h
    TPZSSpStructMatrix.h
    TPZSkylineNSymStructMatrix.h
    TPBSpStructMatrix.h
    TPZParFrontStructMatrix.h
    TPZFrontStructMatrix.h
    pzfstrmatrix.h
    tpzsparseblockdiagonalstructmatrix.h
    #other
    TPZEquationFilter.h
    )
set(headers
    TPBSpStructMatrix.h
    TPZParFrontStructMatrix.h
    TPZSpStructMatrix.h
    pzbdstrmatrix.h
    pzsbstrmatrix.h
    TPZEquationFilter.h
    TPZSSpStructMatrix.h
    TPZStructMatrixBase.h
    pzbstrmatrix.h
    pzskylstrmatrix.h
    pzstrmatrixot.h
    TPZFrontStructMatrix.h
    TPZSkylineNSymStructMatrix.h
    doxstrmatrix.h
    pzfstrmatrix.h
    pzstrmatrix.h
    tpzsparseblockdiagonalstructmatrix.h
    pzstrmatrixflowtbb.h
    TPZStructMatrixTBBFlowUtils.h
   )

set(sources
    TPBSpStructMatrix.cpp
    TPZParFrontStructMatrix.cpp
    TPZSpStructMatrix.cpp
    pzbstrmatrix.cpp
    pzskylstrmatrix.cpp
    pzstrmatrixot.cpp
    TPZEquationFilter.cpp
    TPZSSpStructMatrix.cpp
    TPZStructMatrixBase.cpp
    pzfstrmatrix.cpp
    pzstrmatrix.cpp
    tpzsparseblockdiagonalstructmatrix.cpp
    TPZFrontStructMatrix.cpp
    TPZSkylineNSymStructMatrix.cpp
    pzbdstrmatrix.cpp
    pzsbstrmatrix.cpp
    pzstrmatrixflowtbb.cpp
    TPZStructMatrixTBBFlowUtils.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/StrMatrix)

target_sources(pz PRIVATE ${headers} ${sources})
