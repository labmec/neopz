# @file neopz/Shape/CMakeLists.txt  -- CMake file for the Shape module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Shape>
                           )

set(headers
    TPZShapeData.h
    TPZShapeDisc.h
    pzgenericshape.h
    pzshapeextend.h
    pzshapepiram.h
    pzshapepoint.h
    pzshapequad.h
    pzshapetriang.h
    doxshape.h
    pzshapecube.h
    pzshapelinear.h
    pzshapepiramHdiv.h
    pzshapeprism.h
    pzshapetetra.h
    TPZShapeH1.h
    TPZShapeHDiv.h
    TPZShapeHDivKernel2D.h
    TPZShapeHDivConstant.h
    TPZShapeHDivConstantBound.h
    TPZShapeHCurl.h
    TPZShapeHCurlNoGrads.h
    TPZShapeHDivBound.h
    TPZShapeHDivCollapsed.h
    TPZEnumApproxFamily.h
    )

set(sources
    TPZShapeData.cpp
    TPZShapeDisc.cpp
    pzgenericshape.cpp
    pzshapecube.cpp
    pzshapelinear.cpp
    pzshapepiram.cpp
    pzshapepiramHdiv.cpp
    pzshapeprism.cpp
    pzshapequad.cpp
    pzshapetetra.cpp
    pzshapetriang.cpp
    TPZShapeH1.cpp
    TPZShapeHDiv.cpp
    TPZShapeHDivKernel2D.cpp
    TPZShapeHDivConstant.cpp
    TPZShapeHDivConstantBound.cpp
    TPZShapeHCurl.cpp
    TPZShapeHCurlNoGrads.cpp
    TPZShapeHDivBound.cpp
    TPZShapeHDivCollapsed.cpp
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Shape)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)