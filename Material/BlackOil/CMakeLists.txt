# @file neopz/Material/BlackOil/CMakeLists.txt

set(public_headers
    TPZBlackOil2P3D.h
    )

set(headers
    TPZBlackOil2P3D.h
   )
set(sources
    TPZBlackOil2P3D.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material/BlackOil)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)