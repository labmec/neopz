# @file neopz/Pre/CMakeLists.txt  -- CMake file for the Pre module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Pre>
                           )


set(headers
    MMeshType.h
    TPZExtendGridDimension.h    
    TPZGenSpecialGrid.h
    TPZMHMeshControl.h
    TPZReadGIDGrid.h
    pzhyperplane.h
    TPZAcademicGeoMesh.h
    TPZGMSHReadMesh.h
    TPZGeoMeshBuilder.h
    TPZMHMixedHybridMeshControl.h
    doxpre.h
    pzidentifyrefpattern.h
    pzreadtetgen.h
    TPZAnalyticSolution.h
    TPZGenGrid2D.h
    TPZGmshReader.h
    TPZMHMixedMeshChannelControl.h
    pzbuildmultiphysicsmesh.h
    pzpargrid.h
    tpzhierarquicalgrid.h
    TPZGenGrid3D.h
    TPZHybridizeHDiv.h
    TPZMHMixedMeshControl.h
    pzdatafi.h
    pzreadmesh.h
    )

set(sources
    MMeshType.cpp
    TPZExtendGridDimension.cpp
    TPZGenSpecialGrid.cpp
    TPZMHMeshControl.cpp
    TPZReadGIDGrid.cpp
    pzpargrid.cpp
    TPZAcademicGeoMesh.cpp
    TPZGMSHReadMesh.cpp
    TPZGeoMeshBuilder.cpp
    TPZMHMixedHybridMeshControl.cpp
    pzbuildmultiphysicsmesh.cpp
    TPZAnalyticSolution.cpp
    TPZGenGrid2D.cpp
    TPZGmshReader.cpp
    TPZMHMixedMeshChannelControl.cpp
    pzhyperplane.cpp
    pzreadtetgen.cpp
    TPZGenGrid3D.cpp
    TPZHybridizeHDiv.cpp
    TPZMHMixedMeshControl.cpp
    pzidentifyrefpattern.cpp
    tpzhierarquicalgrid.cpp
    )

if (USING_MKL)
	list(APPEND headers TPZBuildSBFem.h TPZBuildSBFemMultiphysics.h)
	list(APPEND sources TPZBuildSBFem.cpp TPZBuildSBFemMultiphysics.cpp)
endif()

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Pre)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)
