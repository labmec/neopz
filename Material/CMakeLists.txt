# @file neopz/Material/CMakeLists.txt  -- CMake file for the Material module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Material>
                           )

set(public_headers
    TPZMaterial.h
    TPZMaterialT.h
    TPZMatBase.h
    TPZBndCond.h
    TPZBndCondT.h
    TPZBndCondBase.h
    TPZMaterialData.h
    TPZMaterialDataT.h
    TPZMatTypes.h
    #generalinterfaces
    TPZMatError.h
    TPZMatWithMem.h
    TPZMatLoadCases.h
    TPZMatGeneralisedEigenVal.h
    #single space materials
    TPZMatSingleSpace.h
    TPZMatErrorSingleSpace.h
    TPZMatInterfaceSingleSpace.h
    TPZMatTransientSingleSpace.h
    #combined spaces materials
    TPZMatCombinedSpaces.h
    TPZMatErrorCombinedSpaces.h
    TPZMatInterfaceCombinedSpaces.h
    #other
    TPZNullMaterial.h
    TPZLagrangeMultiplier.h
    TPZLagrangeMultiplierCS.h
    )

set(headers
    TPZMaterial.h
    TPZMaterialT.h
    TPZMatBase.h
    TPZBndCond.h
    TPZBndCondT.h
    TPZBndCondBase.h
    TPZMaterialData.h
    TPZMaterialDataT.h
    TPZMatTypes.h
    #generalinterfaces
    TPZMatError.h
    TPZMatWithMem.h
    TPZMatLoadCases.h
    TPZMatGeneralisedEigenVal.h
    #single space materials
    TPZMatSingleSpace.h
    TPZMatErrorSingleSpace.h
    TPZMatInterfaceSingleSpace.h
    TPZMatTransientSingleSpace.h
    #combined spaces materials
    TPZMatCombinedSpaces.h
    TPZMatErrorCombinedSpaces.h
    TPZMatInterfaceCombinedSpaces.h
    #other
    TPZNullMaterial.h
    TPZLagrangeMultiplier.h
    TPZLagrangeMultiplierCS.h
   )
set(sources
    TPZMaterial.cpp
    TPZMaterialT.cpp
    TPZBndCond.cpp
    TPZBndCondT.cpp
    TPZMaterialData.cpp
    TPZMaterialDataT.cpp
    #general interfaces
    TPZMatError.cpp
    TPZMatLoadCases.cpp
    TPZMatGeneralisedEigenVal.cpp
    #single space materials
    TPZMatSingleSpace.cpp
    TPZMatErrorSingleSpace.cpp
    TPZMatInterfaceSingleSpace.cpp
    TPZMatTransientSingleSpace.cpp
    #combined spaces materials
    TPZMatCombinedSpaces.cpp
    TPZMatErrorCombinedSpaces.cpp
    TPZMatInterfaceCombinedSpaces.cpp
    #other
    TPZNullMaterial.cpp
    TPZLagrangeMultiplier.cpp
    TPZLagrangeMultiplierCS.cpp
   )

install(FILES ${public_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material)

target_sources(pz PRIVATE ${headers} ${sources})

#for doxygen
foreach(header ${public_headers})
    list(APPEND PZ_PUBLIC_HEADERS " ${CMAKE_CURRENT_SOURCE_DIR}/${header}")
endforeach()

#materials for conservation laws
add_subdirectory(ConsLaw)
#materials for black oil
add_subdirectory(BlackOil)
#projection materials
add_subdirectory(Projection)
#poisson materials
add_subdirectory(Poisson)

set(PZ_PUBLIC_HEADERS ${PZ_PUBLIC_HEADERS} PARENT_SCOPE)