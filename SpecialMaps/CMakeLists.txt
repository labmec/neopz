# @file neopz/SpecialMaps/CMakeLists.txt  -- CMake file for the SpecialMaps module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/SpecialMaps>
                          )

set(headers
    TPZCylinder.h
    TPZQuadTorus.h
    TPZTriangleTorus.h
    convtest.h
    tpzblendnaca.h
    tpzellipse3d.h
    tpzmathtools.h
    tpzquadraticline.h
    tpzquadraticpyramid.h
    tpzquadratictetra.h
    TPZQuadSphere.h
    TPZTriangleSphere.h
    TPZWavyLine.h
    tpzarc3d.h
    tpzchangeel.h
    tpzgeomid.h
    tpzquadraticcube.h
    tpzquadraticprism.h
    tpzquadraticquad.h
    tpzquadratictrig.h
    )

set(sources
    TPZCylinder.cpp
    TPZQuadTorus.cpp
    TPZTriangleTorus.cpp
    convtest.cpp
    tpzblendnaca.cpp
    tpzellipse3d.cpp
    tpzquadraticcube.cpp
    tpzquadraticprism.cpp
    tpzquadraticquad.cpp
    tpzquadratictrig.cpp
    TPZQuadSphere.cpp
    TPZTriangleSphere.cpp
    TPZWavyLine.cpp
    tpzarc3d.cpp
    tpzchangeel.cpp
    tpzmathtools.cpp
    tpzquadraticline.cpp
    tpzquadraticpyramid.cpp
    tpzquadratictetra.cpp
    )

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/SpecialMaps)

target_sources(pz PRIVATE ${headers} ${sources})
