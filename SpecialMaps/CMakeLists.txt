# @file neopz/SpecialMaps/CMakeLists.txt  -- CMake file for the SpecialMaps module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SpecialMaps>
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

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SpecialMaps)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)