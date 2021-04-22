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
    #full
    pzfstrmatrix.h
    pzsfstrmatrix.h
    #band
    pzbstrmatrix.h
    pzsbstrmatrix.h
    #skyline
    TPZSkylineNSymStructMatrix.h
    pzskylstrmatrix.h
    #sparse
    TPZSpStructMatrix.h
    TPZSSpStructMatrix.h
    #front
    TPZFrontStructMatrix.h
    TPZParFrontStructMatrix.h
    #block
    pzbdstrmatrix.h
    TPZBSpStructMatrix.h
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
    #full
    pzfstrmatrix.h
    pzsfstrmatrix.h
    #band
    pzbstrmatrix.h
    pzsbstrmatrix.h
    #skyline
    TPZSkylineNSymStructMatrix.h
    pzskylstrmatrix.h
    #sparse
    TPZSpStructMatrix.h
    TPZSSpStructMatrix.h
    #front
    TPZFrontStructMatrix.h
    TPZParFrontStructMatrix.h
    #block
    pzbdstrmatrix.h
    TPZBSpStructMatrix.h
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
    #full
    pzfstrmatrix.cpp
    pzsfstrmatrix.cpp
    #band
    pzbstrmatrix.cpp
    pzsbstrmatrix.cpp
    #skyline
    TPZSkylineNSymStructMatrix.cpp
    pzskylstrmatrix.cpp
    #sparse
    TPZSpStructMatrix.cpp
    TPZSSpStructMatrix.cpp
    #front
    TPZFrontStructMatrix.cpp
    TPZParFrontStructMatrix.cpp
    #block
    pzbdstrmatrix.cpp
    TPZBSpStructMatrix.cpp
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