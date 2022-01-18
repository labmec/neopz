# @file neopz/Material/Projection/CMakeLists.txt

set(public_headers
    TPZL2Projection.h
    TPZHDivProjection.h
    TPZHCurlProjection.h
    TPZL2ProjectionCS.h
    )

set(headers
    TPZL2Projection.h
    TPZHDivProjection.h
    TPZHCurlProjection.h
    TPZL2ProjectionCS.h
   )
set(sources
    TPZL2Projection.cpp
    TPZHDivProjection.cpp
    TPZHCurlProjection.cpp
    TPZL2ProjectionCS.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material/Projection)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)
