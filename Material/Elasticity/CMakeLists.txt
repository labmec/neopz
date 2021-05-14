# @file neopz/Material/Elasticity/CMakeLists.txt

set(public_headers
    TPZElasticity3D.h
    TPZElasticity2D.h
    )

set(headers
    TPZElasticity3D.h
    TPZElasticity2D.h
   )
set(sources
    TPZElasticity3D.cpp
    TPZElasticity2D.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material/Elasticity)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)