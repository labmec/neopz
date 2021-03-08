# @file neopz/Frontal/CMakeLists.txt  -- CMake file for the Frontal module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/Frontal>
                          )

set(headers
    TPZFileEqnStorage.h
    TPZFront.h
    TPZFrontMatrix.h
    TPZFrontNonSym.h
    TPZFrontSym.h
    TPZParFrontMatrix.h
    TPZStackEqnStorage.h
    tpzeqnarray.h
    )

set(sources
    TPZFileEqnStorage.cpp
    TPZFront.cpp
    TPZFrontMatrix.cpp
    TPZFrontNonSym.cpp
    TPZFrontSym.cpp
    TPZParFrontMatrix.cpp
    TPZStackEqnStorage.cpp
    pzvecfrontal.cpp
    tpzeqnarray.cpp
    )

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Frontal)

target_sources(pz PRIVATE ${headers} ${sources})
