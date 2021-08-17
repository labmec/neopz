# @file neopz/Util/CMakeLists.txt  -- CMake file for the Util module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Util>
                           )
set(public_headers
  #containers
  pzvec.h
  pzmanvector.h
  TPZChunkVector.h
  pzadmchunk.h
  pzstack.h
  #parallel
  TPZParallelUtils.h
  TPZSemaphore.h
  #other
  TPZTimer.h
  TPZSimpleTimer.h
  tpzpermutation.h
  pzaxestools.h
  pzlog.h
  pzextractval.h
  tpzautopointer.h  
  pzfunction.h
  TPZExactFunction.h
  pzvec_extras.h
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
    TPZSimpleTimer.h
    pzadmchunk.h
    pzcheckconsistency.h
    pzgradient.h
    pzmanvector.h
    pzpix.h
    pzstack.h
    pzvec_extras.h
    TPZExactFunction.h
    )

set(sources
    TPZAssert.cpp
    TPZParallelUtils.cpp
    pzcheckconsistency.cpp
    pzline.cpp
    pznuma.cpp
    pzpolynomial.cpp
    TPZFunctionTranslator.cpp
    TPZSemaphore.cpp
    pzfunction.cpp
    pzlog.cpp
    pznumeric.cpp
    tpzpermutation.cpp
    TPZHWTools.cpp
    TPZTimer.cpp
    TPZSimpleTimer.cpp
    pzgradient.cpp
    pzplane.cpp
    pzstring.cpp
    tpzautopointer.cpp
    pzvecgen.cpp
    pzmanvectorgen.cpp
    pzstackgen.cpp
    )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Util)

if(USING_LOG4CXX)
    set(inputPath ${CMAKE_CURRENT_SOURCE_DIR}/log4cxx.cfg)
    set(outputPath ${CMAKE_CURRENT_BINARY_DIR}/log4cxx.cfg)

    add_custom_command(
        OUTPUT ${outputPath}
        DEPENDS ${inputPath}
        COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/cmake/log4cxx_config_backup_gen.cmake
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${inputPath} ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_custom_target(log4cxx-config DEPENDS ${outputPath})
    add_dependencies(pz log4cxx-config)

    install(FILES ${outputPath} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Util)
endif()

target_sources(pz PRIVATE ${headers} ${sources})
add_subdirectory(Hash)
#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)
