# @file neopz/Material/Poisson/CMakeLists.txt

set(public_headers
    TPZWgma.h
    TPZPlanarWgScatt.h
    TPZPlanarWgScattSrc.h
    TPZPeriodicWgma.h
    TPZMatPML.h
    TPZMatPML_impl.h
    )

set(headers
    TPZWgma.h
    TPZPlanarWgScatt.h
    TPZPlanarWgScattSrc.h
    TPZPeriodicWgma.h
    TPZMatPML.h
    TPZMatPML_impl.h
   )
set(sources
    TPZWgma.cpp
    TPZPlanarWgScatt.cpp
    TPZPlanarWgScattSrc.cpp
    TPZPeriodicWgma.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material/Electromagnetics)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)