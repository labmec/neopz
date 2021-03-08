# @file neopz/Util/CMakeLists.txt  -- CMake file for the Util module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/Util>
                          )

set(headers
    TPZAssert.h
    TPZLimitedPriorityQueue.h
    TPZTimer.h
    pzadmchunkthreadsafe.h
    pzfunction.h
    pzmanvector.h
    pzplane.h
    pzvec.h
    tpzpermutation.h
    TPZChunkVector.h
    TPZParallelUtils.h
    checkconv.h
    pzaxestools.h
    pzgradient.h
    pznuma.h
    pzpolynomial.h
    pzvec_extras.h
    TPZFunctionTranslator.h
    TPZPriorityQueue.h
    doxutil.h
    pzcheckconsistency.h
    pzline.h
    pznumeric.h
    pzstack.h
    tpzautopointer.h
    )

set(sources
    TPZAssert.cpp
    TPZParallelUtils.cpp
    pzcheckconsistency.cpp
    pzline.cpp
    pznuma.cpp
    pzpolynomial.cpp
    pzvecgen.cpp
    tpzpagemigrationmanager.cpp
    TPZFunctionTranslator.cpp
    TPZSemaphore.cpp
    pzfunction.cpp
    pzlog.cpp
    pznumeric.cpp
    pzstackgen.cpp
    pzvecmesh.cpp
    tpzpermutation.cpp
    TPZHWTools.cpp
    TPZTimer.cpp
    pzgradient.cpp
    pzmanvectorgen.cpp
    pzplane.cpp
    pzstring.cpp
    tpzautopointer.cpp
    )

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Util)

target_sources(pz PRIVATE ${headers} ${sources})

add_subdirectory(Hash)
