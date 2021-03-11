
target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Random>
                           )

set(headers
    TPZConstrainedNormalRandom.h
    TPZConstrainedRandom.h
    TPZNormalRandom.h
    TPZRandom.h
    TPZUniformRandom.h
    )

set(sources TPZUniformRandom.cpp)

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Random)

target_sources(pz PRIVATE ${headers} ${sources})
