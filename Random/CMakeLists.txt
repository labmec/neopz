
target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${NEOPZ_INCLUDE}/Random>
                          )

set(headers
    TPZConstrainedNormalRandom.h
    TPZConstrainedRandom.h
    TPZNormalRandom.h
    TPZRandom.h
    TPZUniformRandom.h
    )

set(sources TPZUniformRandom.cpp)

install(FILES ${headers} DESTINATION ${NEOPZ_INCLUDE}/Random)

target_sources(pz PRIVATE ${headers} ${sources})
