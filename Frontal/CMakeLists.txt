# @file neopz/Frontal/CMakeLists.txt  -- CMake file for the Frontal module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Frontal>
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

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Frontal)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)