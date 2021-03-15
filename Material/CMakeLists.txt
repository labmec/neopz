# @file neopz/Material/CMakeLists.txt  -- CMake file for the Material module

target_include_directories(pz PUBLIC 
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Material>
                           )

set(headers
    TPZBndCondWithMem.h
    TPZLagrangeMultiplier.h
    TPZMatLaplacian.h
    TPZNullMaterial.h
    TPZReynoldsFlow.h
    pzbndcond.h
    pzl2projection.h
    pzmaterialdata.h
    tpzoutofrange.h
    TPZBndCondWithMem_impl.h
    TPZMatHCurlProjection.h
    TPZMaterial.h
    TPZNullMaterialTranslator.h
    TPZVecL2.h
    pzconslaw.h
    pzmat1dlin.h
    pzpoisson3d.h
    TPZMatHelmholtz.h
    TPZMaterialTranslator.h
    TPZPorousElasticity.h
    doxmaterial.h
    pzmat2dlin.h
    pztransientmat.h
    )

set(sources
    TPZBndCondWithMem.cpp
    TPZMatHCurlProjection.cpp
    TPZMaterial.cpp
    TPZNullMaterialTranslator.cpp
    TPZVecL2.cpp
    pzmat2dlin.cpp
    pztransientmat.cpp
    TPZMatHelmholtz.cpp
    TPZMaterialTranslator.cpp
    TPZPorousElasticity.cpp
    pzbndcond.cpp
    pzl2projection.cpp
    pzmaterialdata.cpp
    tpzoutofrange.cpp
    TPZLagrangeMultiplier.cpp
    TPZMatLaplacian.cpp
    TPZNullMaterial.cpp
    TPZReynoldsFlow.cpp
    pzconslaw.cpp
    pzmat1dlin.cpp
    pzpoisson3d.cpp
    )

install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material)

target_sources(pz PRIVATE ${headers} ${sources})
if (BUILD_PLASTICITY_MATERIALS)
	if ((NOT STATE_TYPE STREQUAL "double") OR (NOT REAL_TYPE STREQUAL "double"))
		message (FATAL_ERROR
      "\nConflicting options for BUILD_PLASTICITY_MATERIALS\n"
      "Set REAL_TYPE=double and STATE_TYPE=double to build them\n"
      "Or set BUILD_PLASTICITY_MATERIALS=OFF")
  endif()
endif()
if(BUILD_COMPLEX_PROJECTS)
  add_subdirectory(Complex)
else()
  add_subdirectory(REAL)
endif()
