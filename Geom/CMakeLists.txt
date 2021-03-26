# @file neopz/Geom/CMakeLists.txt  -- CMake file for the Geom module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Geom>
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

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Geom)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)