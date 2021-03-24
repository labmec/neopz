# @file neopz/Topology/CMakeLists.txt  -- CMake file for the Topology module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Topology>
                           )

set(headers
    PrismExtend.h
    TPZTopologyUtils.h
    doxtopology.h
    tpzcube.h
    tpzline.h
    tpzpoint.h
    tpzprism.h
    tpzpyramid.h
    tpzquadrilateral.h
    tpztetrahedron.h
    tpztriangle.h
    )

set(sources
    PrismExtend.cpp
    TPZTopologyUtils.cpp
    tpzcube.cpp
    tpzline.cpp
    tpzpoint.cpp
    tpzprism.cpp
    tpzpyramid.cpp
    tpzquadrilateral.cpp
    tpztetrahedron.cpp
    tpztriangle.cpp
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Topology)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)