# @file neopz/StrMatrix/CMakeLists.txt  -- CMake file for the StrMatrix module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/StrMatrix>
                           )
#@TODO: check if these public headers are exaggerated
set(public_headers
    TPZStrMatParInterface.h
    TPZStructMatrix.h
    TPZStructMatrixT.h
    #parallel scheme
    pzstrmatrixor.h
    pzstrmatrixot.h
    pzstrmatrixflowtbb.h
    #matrix storage format
    pzbstrmatrix.h
    TPZSpStructMatrix.h
    pzbdstrmatrix.h
    pzsbstrmatrix.h
    pzskylstrmatrix.h
    TPZSSpStructMatrix.h
    TPZSkylineNSymStructMatrix.h
    TPZBSpStructMatrix.h
    TPZParFrontStructMatrix.h
    TPZFrontStructMatrix.h
    pzfstrmatrix.h
    # tpzsparseblockdiagonalstructmatrix.h # needsreview
    #other
    TPZEquationFilter.h
   )
set(headers
    TPZStrMatParInterface.h
    TPZStructMatrix.h
    TPZStructMatrixT.h
    #parallel scheme
    pzstrmatrixor.h
    pzstrmatrixot.h
    pzstrmatrixflowtbb.h
    #matrix storage format
    pzbstrmatrix.h
    TPZSpStructMatrix.h
    pzbdstrmatrix.h
    pzsbstrmatrix.h
    pzskylstrmatrix.h
    TPZSSpStructMatrix.h
    TPZSkylineNSymStructMatrix.h
    TPZBSpStructMatrix.h
    TPZParFrontStructMatrix.h
    TPZFrontStructMatrix.h
    pzfstrmatrix.h
    tpzsparseblockdiagonalstructmatrix.h
    #other
    TPZEquationFilter.h
    TPZStructMatrixTBBFlowUtils.h
   )

set(sources
    TPZStrMatParInterface.cpp
    TPZStructMatrix.cpp
    TPZStructMatrixT.cpp
    #parallel schemes
    pzstrmatrixor.cpp
    pzstrmatrixot.cpp
    pzstrmatrixflowtbb.cpp
    #matrix storage format
    pzbstrmatrix.cpp
    TPZSpStructMatrix.cpp
    pzbdstrmatrix.cpp
    pzsbstrmatrix.cpp
    pzskylstrmatrix.cpp
    TPZSSpStructMatrix.cpp
    TPZSkylineNSymStructMatrix.cpp
    TPZBSpStructMatrix.cpp
    TPZParFrontStructMatrix.cpp
    TPZFrontStructMatrix.cpp
    pzfstrmatrix.cpp
    tpzsparseblockdiagonalstructmatrix.cpp
    #other
    TPZEquationFilter.cpp
    TPZStructMatrixTBBFlowUtils.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/StrMatrix)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)