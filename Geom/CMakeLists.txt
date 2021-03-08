# @file neopz/Geom/CMakeLists.txt  -- CMake file for the Geom module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/Geom>
                          )

set(headers
    TPZGeoCube.h
    TPZGeoLinear.h
    pzgeom_utility.h
    pzgeoprism.h
    pzgeoquad.h
    pzgeotriangle.h
    pznoderep.h.h
    tpzgeoblend.h
    TPZGeoExtend.h
    doxgeometry.h
    pzgeopoint.h
    pzgeopyramid.h
    pzgeotetrahedra.h
    pznoderep.h
    pzshapeextend.h
    )

set(sources
    TPZGeoCube.cpp
    TPZGeoLinear.cpp
    pzgeom_utility.cpp
    pzgeopoint.cpp
    pzgeoprism.cpp
    pzgeopyramid.cpp
    pzgeoquad.cpp
    pzgeotetrahedra.cpp
    pzgeotriangle.cpp
    pznoderep.cpp
    tpzgeoblend.cpp
    )

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Geom)

target_sources(pz PRIVATE ${headers} ${sources})
