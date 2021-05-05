# @file neopz/Material/ConsLaw/CMakeLists.txt

set(public_headers
    TPZConsLaw.h
    )

set(headers
    TPZConsLaw.h
   )
set(sources
    TPZConsLaw.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material/ConsLaw)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)