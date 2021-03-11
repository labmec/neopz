# @file neopz/StrMatrix/CMakeLists.txt  -- CMake file for the StrMatrix module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/StrMatrix>
                           )

set(headers
    TPBSpStructMatrix.h
    TPZParFrontStructMatrix.h
    TPZSpStructMatrix.h
    pzbdstrmatrix.h
    pzsbstrmatrix.h
    pzstrmatrixflowtbb.h
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
   )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/StrMatrix)

target_sources(pz PRIVATE ${headers} ${sources})
