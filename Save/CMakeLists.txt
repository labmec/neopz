# @file neopz/Save/CMakeLists.txt  -- CMake file for the Save module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Save>
                           )
#@TODO: define public interface
set(headers
    TPZBFileStream.h
    TPZChunkTranslator.h
    TPZContBufferedStream.h
    TPZGeneralFStream.h
    TPZRestoredInstance.h
    TPZStream.h
    pzmd5stream.h
    TPZChunkInTranslation.h
    TPZCircBufferedStream.h
    TPZFileStream.h
    TPZPersistenceManager.h
    TPZSavable.h
    doxsave.h
    )

set(sources
    TPZBFileStream.cpp
    TPZChunkTranslator.cpp
    TPZContBufferedStream.cpp
    TPZPersistenceManager.cpp
    TPZSavable.cpp
    pzmd5stream.cpp
    TPZChunkInTranslation.cpp
    TPZCircBufferedStream.cpp
    TPZFileStream.cpp
    TPZRestoredInstance.cpp
    TPZStream.cpp
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Save)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)