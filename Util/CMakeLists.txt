# @file neopz/Util/CMakeLists.txt  -- CMake file for the Util module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Util>
                           )

set(headers
    TPZAssert.h
    TPZHWTools.h
    TPZPriorityQueue.h
    checkconv.h
    pzadmchunkthreadsafe.h
    pzextractval.h
    pzline.h
    pznuma.h
    pzplane.h
    pzstring.h
    tpzautopointer.h
    TPZChunkVector.h
    TPZLimitedPriorityQueue.h
    TPZSemaphore.h
    doxutil.h
    pzaxestools.h
    pzfunction.h
    pzlog.h
    pznumeric.h
    pzpolynomial.h
    pzvec.h
    tpzpermutation.h
    TPZFunctionTranslator.h
    TPZParallelUtils.h
    TPZTimer.h
    pzadmchunk.h
    pzcheckconsistency.h
    pzgradient.h
    pzmanvector.h
    pzpix.h
    pzstack.h
    pzvec_extras.h
    )

set(sources
    TPZAssert.cpp
    TPZParallelUtils.cpp
    pzcheckconsistency.cpp
    pzline.cpp
    pznuma.cpp
    pzpolynomial.cpp
    pzvecgen.cpp
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

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Util)

target_sources(pz PRIVATE ${headers} ${sources})

add_subdirectory(Hash)
