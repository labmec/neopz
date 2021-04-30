# @file neopz/SubStruct/CMakeLists.txt  -- CMake file for the SubStruct module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SubStruct>
                           )
set(public_headers
    doxsubstrmatrix.h
    tpzdohrassemblelist.h
    tpzdohrmatrix.h
    tpzdohrsubstruct.h
    tpzgensubstruct.h
    # tpzpairstructmatrix.h #needsreview
    tpzdohrassembly.h
    tpzdohrprecond.h
    tpzdohrsubstructCondense.h
    tpzmatredstructmatrix.h
    pzdohrstructmatrix.h
    )
set(headers
    TPZTimeTemp.h
    doxsubstrmatrix.h
    tpzdohrassemblelist.h
    tpzdohrmatrix.h
    tpzdohrsubstruct.h
    tpzgensubstruct.h
    tpzpairstructmatrix.h
    tpzdohrassembly.h
    tpzdohrprecond.h
    tpzdohrsubstructCondense.h
    tpzmatredstructmatrix.h
    tpzparallelenviroment.h
    pzdohrstructmatrix.h
    )

set(sources
    TPZTimeTemp.cpp
    tpzdohrassembly.cpp
    tpzdohrprecond.cpp
    tpzdohrsubstructCondense.cpp
    tpzmatredstructmatrix.cpp
    tpzparallelenviroment.cpp
    tpzdohrassemblelist.cpp
    tpzdohrmatrix.cpp
    tpzdohrsubstruct.cpp
    tpzgensubstruct.cpp
    tpzpairstructmatrix.cpp
    pzdohrstructmatrix.cpp
    )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SubStruct)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)
