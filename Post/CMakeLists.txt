# @file neopz/Post/CMakeLists.txt  -- CMake file for the Post module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Post>
                          )

set(headers
    TPZDrawStyle.h
    TPZVTKGeoMesh.h
    pzdxmesh.h
    pzgraphel1d.h
    pzgraphelq2dd.h
    pzgraphnode.h
    pzv3dmesh.h
    tpzgraphelprismmapped.h
    tpzgraphelt3d.h
    TPZMeshSolution.h
    doxpost.h
    pzgraphel1dd.h
    pzgraphelq3dd.h
    pzmvmesh.h
    pztrigraph.h
    pzvisualmatrix.h
    tpzgraphelpyramidmapped.h
    TPZProjectEllipse.h
    pzgraphel.h
    pzgraphelq2d.h
    pzgraphmesh.h
    pztrigraphd.h
    pzvtkmesh.h
    tpzgraphelt2dmapped.h
    pzgradientreconstruction.h
    #post proc
    pzpostprocanalysis.h
    pzcompelpostproc.h
    pzpostprocmat.h    
    )

set(sources
    TPZMeshSolution.cpp
    pzdxmesh.cpp
    pzgraphel1d.cpp
    pzgraphelq2dd.cpp
    pzgraphnode.cpp
    pzv3dmesh.cpp
    tpzgraphelprismmapped.cpp
    tpzgraphelt3d.cpp
    TPZProjectEllipse.cpp
    pzgraphel1dd.cpp
    pzgraphelq3dd.cpp
    pzmvmesh.cpp
    pztrigraph.cpp
    pzvisualmatrix.cpp
    tpzgraphelpyramidmapped.cpp
    TPZVTKGeoMesh.cpp
    pzgraphel.cpp
    pzgraphelq2d.cpp
    pzgraphmesh.cpp
    pztrigraphd.cpp
    pzvtkmesh.cpp
    tpzgraphelt2dmapped.cpp
    pzgradientreconstruction.cpp
    pzpostprocanalysis.cpp
    pzpostprocmat.cpp
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Post)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)