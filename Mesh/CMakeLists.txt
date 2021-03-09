# @file neopz/Mesh/CMakeLists.txt  -- CMake file for the Mesh module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/Mesh>
                          )


set(headers
  doxmesh.h
  pzconnect.h
  pzgeoelbc.h  
  TPZCompMeshTools.h
  TPZGeoMeshTools.h
  pzcheckgeom.h
  pzcreateapproxspace.h
  pzgeoel.h
  pzintel.h  
  TPZGeoElement.h
  TPZInterfaceEl.h
  pzcheckmesh.h
  pzelchdivbound2.h
  pzgeoelrefless.h
  pzinterpolationspace.h
  TPZCompElDisc.h
  TPZGeoElement.h.h
  TPZMultiphysicsCompMesh.h
  pzcheckrestraint.h
  pzelchdiv.h
  pzgeoelrefless.h.h
  pzmultiphysicscompel.h
  TPZCompElHCurlFull.h
  tpzgeoelmapped.h
  TPZMultiphysicsInterfaceEl.h
  pzcmesh.h
  pzelctemp.h
  pzgeoelside.h
  pzmultiphysicselement.h
  TPZCompElHCurl.h
  tpzgeoelrefpattern.h
  TPZOneShapeRestraint.h
  pzcompel.h
  pzelementgroup.h
  pzgmesh.h
  pzreducedspace.h
  TPZCompElHDivCollapsed.h
  tpzgeoelrefpattern.h.h
  TPZSBFemElementGroup.h
  pzcompelwithmem.h
  pzelmat.h
  pzgnode.h
  pzreferredcompel.h
  TPZCompElLagrange.h
  TPZGeoElSideAncestors.h
  TPZSBFemVolume.h
  pzcondensedcompel.h    
  pzsubcmesh.h
  tpzcompmeshreferred.h
  TPZGeoElSidePartition.h
)

set(sources
  pzcheckgeom.cpp
  pzelchdivbound2.cpp
  pzgeoelside.cpp
  pzreducedspace.cpp
  TPZCompElHCurlFull.cpp
  TPZGeoElSidePartition.cpp
  pzcheckmesh.cpp
  pzelchdiv.cpp
  pzgmesh.cpp
  pzreferredcompel.cpp
  TPZCompElHDivCollapsed.cpp
  TPZGeoMeshTools.cpp
  pzcheckrestraint.cpp
  pzelctemp.cpp
  pzgnode.cpp
  pzstackmesh.cpp
  TPZCompElLagrange.cpp
  TPZInterfaceEl.cpp
  pzcmesh.cpp
  pzelementgroup.cpp  
  pzsubcmesh.cpp
  tpzcompmeshreferred.cpp
  TPZMultiphysicsCompMesh.cpp
  pzcompel.cpp
  pzelmat.cpp  
  pzvecmesh.cpp
  TPZCompMeshTools.cpp
  TPZMultiphysicsInterfaceEl.cpp
  pzcompelwithmem.cpp  
  pzintel.cpp  
  TPZGeoElement.cpp
  TPZOneShapeRestraint.cpp
  pzcondensedcompel.cpp
  pzgeoelbc.cpp
  pzinterpolationspace.cpp
  TPZAgglomerateMesh.cpp
  tpzgeoelmapped.cpp
  TPZSBFemElementGroup.cpp
  pzconnect.cpp
  pzgeoel.cpp
  pzmultiphysicscompel.cpp
  TPZCompElDisc.cpp
  tpzgeoelrefpattern.cpp
  TPZSBFemVolume.cpp
  pzcreateapproxspace.cpp
  pzgeoelrefless.cpp
  pzmultiphysicselement.cpp
  TPZCompElHCurl.cpp
  TPZGeoElSideAncestors.cpp
)
if(NOT BUILD_COMPLEX_PROJETS)
  list(APPEND headers pzflowcmesh.h pzhdivpressure.h 
    pzhdivpressurebound.h tpzagglomeratemesh.h TPZAgglomerateEl.h)
  list(APPEND sources pzflowcmesh.cpp pzhdivpressure.cpp
    pzhdivpressurebound.cpp TPZAgglomerateEl.cpp)
  
  if(USING_MKL)
    list(APPEND headers TPZSBFemElementGroup.h TPZSBFemVolume.h)
    list(APPEND sources TPZSBFemElementGroup.cpp TPZSBFemVolume.cpp)
  endif()
endif()
  
install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Mesh)

target_sources(pz PRIVATE ${headers} ${sources})

