# @file neopz/SubStruct/CMakeLists.txt  -- CMake file for the SubStruct module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/SubStruct>
                          )

set(headers
    TPZTimeTemp.h
    doxsubstrmatrix.h
    tpzdohrassemblelist.h
    tpzdohrmatrix.h
    tpzdohrsubstruct.h
    tpzgensubstruct.h
    tpzpairstructmatrix.h
    TPZfTime.h
    pzdohrstructmatrix.h
    tpzdohrassembly.h
    tpzdohrprecond.h
    tpzdohrsubstructCondense.h
    tpzmatredstructmatrix.h
    tpzparallelenviroment.h
    )

set(sources
    TPZTimeTemp.cpp
    pzdohrstructmatrix.cpp
    tpzdohrassembly.cpp
    tpzdohrprecond.cpp
    tpzdohrsubstructCondense.cpp
    tpzmatredstructmatrix.cpp
    tpzparallelenviroment.cpp
    TPZfTime.cpp
    tpzdohrassemblelist.cpp
    tpzdohrmatrix.cpp
    tpzdohrsubstruct.cpp
    tpzgensubstruct.cpp
    tpzpairstructmatrix.cpp
    )

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/SubStruct)

target_sources(pz PRIVATE ${headers} ${sources})
