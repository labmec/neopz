# @file neopz/Shape/CMakeLists.txt  -- CMake file for the Shape module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/Shape>
                          )

set(headers
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
    )

set(sources
    TPZShapeDisc.cpp
    pzshapecube.cpp
    pzshapelinear.cpp
    pzshapepiram.cpp
    pzshapepiramHdiv.cpp
    pzshapeprism.cpp
    pzshapequad.cpp
    pzshapetetra.cpp
    pzshapetriang.cpp
    )

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Shape)

target_sources(pz PRIVATE ${headers} ${sources})
