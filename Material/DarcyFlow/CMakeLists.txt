# @file neopz/Material/DarcyFlow/CMakeLists.txt

set(public_headers
    TPZDarcyFlowInterface.h
    TPZDarcyFlow.h
    TPZMixedDarcyFlow.h
    )

set(headers
    TPZDarcyFlowInterface.h
    TPZDarcyFlow.h
    TPZMixedDarcyFlow.h
    )
set(sources
    TPZDarcyFlowInterface.cpp
    TPZDarcyFlow.cpp
    TPZMixedDarcyFlow.cpp
    )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material/DarcyFlow)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)
