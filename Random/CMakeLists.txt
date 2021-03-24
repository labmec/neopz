
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

#for doxygen
foreach(header ${headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()
set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)