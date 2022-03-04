# @file neopz/Material/Poisson/CMakeLists.txt

set(public_headers
    TPZWaveguideModalAnalysis.h
    TPZWaveguideModalAnalysisPML.h
    TPZPlanarWGScattering.h
    TPZPlanarWGScatteringPML.h
    )

set(headers
    TPZWaveguideModalAnalysis.h
    TPZWaveguideModalAnalysisPML.h
    TPZPlanarWGScattering.h
    TPZPlanarWGScatteringPML.h
   )
set(sources
    TPZWaveguideModalAnalysis.cpp
    TPZWaveguideModalAnalysisPML.cpp
    TPZPlanarWGScattering.cpp
    TPZPlanarWGScatteringPML.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material/Electromagnetics)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)