# @file neopz/Material/REAL/Plasticity/CMakeLists.txt  -- CMake file for Plasticity Materials

target_include_directories(pz PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
                           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/Material/REAL/Plasticity>
                           )
set(headers
    TPZDruckerPrager.h
    TPZElasticCriterion.h
    TPZElasticResponse.h
    TPZElastoPlasticMem.h
    TPZLadeKim.h
    TPZLadeKimThermoForceA.h
    TPZLadeNelsonElasticResponse.h
    TPZMatElastoPlastic.h
    TPZMatElastoPlastic2D.h
    TPZMatElastoPlastic2D_impl.h
    TPZMatElastoPlastic_impl.h
    TPZMatPorous.h
    TPZMatPorous_impl.h
    TPZMatTemporal.h
    TPZMohrCoulomb.h
    TPZMohrCoulombNeto.h
    TPZMohrCoulombPV.h
    TPZPlasticBase.h
    TPZPlasticCriterion.h
    TPZPlasticIntegrMem.h
    TPZPlasticState.h
    TPZPlasticStep.h
    TPZPlasticStepID.h
    TPZPlasticStepPV.h
    TPZPlasticityTest.h
    TPZPorousElasticCriterion.h
    TPZPorousElasticResponse.h
    TPZPorousElastoPlasticMem.h
    TPZSandlerDimaggio.h
    TPZSandlerDimaggioThermoForceA.h
    TPZSandlerExtended.h
    TPZThermoForceA.h
    TPZVonMises.h
    TPZWillamWarnke.h
    TPZYCCamClayPV.h
    TPZYCDruckerPrager.h
    TPZYCDruckerPragerPV.h
    TPZYCLadeKim.h
    TPZYCModifiedMohrCoulomb.h
    TPZYCMohrCoulomb.h
    TPZYCMohrCoulombPV.h
    TPZYCRankine.h
    TPZYCSandlerDimaggio.h
    TPZYCSandlerDimaggioL.h
    TPZYCSandlerDimaggioL2.h
    TPZYCTresca.h
    TPZYCVonMises.h
    TPZYCWillamWarnke.h
    pzelastoplasticanalysis.h
    pzplasticdiagnostic.h
    pzporoanalysis.h
    pzporoelastoplasticmem.h
    tpzyctrescaregularized.h
    tpzycvonmises_interface.h
    tpzycvonmisescombtresca.h
    )

set(sources
    TPZDruckerPrager.cpp
    TPZElasticCriterion.cpp
    TPZElasticResponse.cpp
    TPZElastoPlasticMem.cpp
    TPZLadeKim.cpp
    TPZLadeKimThermoForceA.cpp
    TPZLadeNelsonElasticResponse.cpp
    TPZMatElastoPlastic.cpp
    TPZMatElastoPlastic2D.cpp
    TPZMatPorous.cpp
    TPZMatTemporal.cpp
    TPZMohrCoulomb.cpp
    TPZPlasticBase.cpp
    TPZPlasticState.cpp
    TPZPlasticStep.cpp
    TPZPlasticStepPV.cpp
    TPZPorousElasticCriterion.cpp
    TPZPorousElasticResponse.cpp
    TPZPorousElastoPlasticMem.cpp
    TPZSandlerDimaggioThermoForceA.cpp
    TPZSandlerExtended.cpp
    TPZThermoForceA.cpp
    TPZVonMises.cpp
    TPZWillamWarnke.cpp
    TPZYCCamClayPV.cpp
    TPZYCDruckerPrager.cpp
    TPZYCDruckerPragerPV.cpp
    TPZYCLadeKim.cpp
    TPZYCModifiedMohrCoulomb.cpp
    TPZYCMohrCoulomb.cpp
    TPZYCMohrCoulombPV.cpp
    TPZYCSandlerDimaggio.cpp
    TPZYCSandlerDimaggioL.cpp
    TPZYCSandlerDimaggioL2.cpp
    TPZYCTresca.cpp
    TPZYCTrescaRegularized.cpp
    TPZYCVonMises.cpp
    TPZYCWillamWarnke.cpp
    pzelastoplasticanalysis.cpp
    pzplasticdiagnostic.cpp
    pzporoanalysis.cpp
    pzporoelastoplasticmem.cpp
    tpzycvonmisescombtresca.cpp
    )

set(PZ_FIX_PLASTIC_TRANSLATORS OFF)
if(PZ_FIX_PLASTIC_TRANSLATORS)
    add_compile_definition(pz PRIVATE FIX_PLASTIC_TRANSLATORS)
    list(APPEND headers
        TPZElasticCriterionTranslator.h
        TPZDruckerPragerTranslator.h
        TPZElasticResponseTranslator.h
        TPZElastoPlasticMemTranslator.h
        TPZMatElastoPlastic2DTranslator.h
        TPZMatElastoPlasticTranslator.h
        TPZPlasticStateTranslator.h
        TPZPlasticStepPVTranslator.h
        TPZPlasticStepTranslator.h
        TPZPorousElasticCriterionTranslator.h
        TPZPorousElasticResponseTranslator.h
        TPZSandlerDimaggioThermoForceATranslator.h
        TPZSandlerDimaggioTranslator.h
        TPZSandlerExtendedTranslator.h    
        TPZYCCamClayPVTranslator.h    
        TPZYCDruckerPragerPVTranslator.h    
        TPZYCMohrCoulombPVTranslator.h    
        TPZYCSandlerDimaggioL2Translator.h
        TPZYCSandlerDimaggioLTranslator.h
        TPZYCSandlerDimaggioTranslator.h
        )

    list(APPEND sources
        TPZElasticCriterionTranslator.cpp
        TPZDruckerPragerTranslator.cpp
        TPZElasticResponseTranslator.cpp
        TPZElastoPlasticMemTranslator.cpp
        TPZMatElastoPlastic2DTranslator.cpp
        TPZMatElastoPlasticTranslator.cpp
        TPZPlasticStateTranslator.cpp
        TPZPlasticStepPVTranslator.cpp
        TPZPlasticStepTranslator.cpp
        TPZPorousElasticCriterionTranslator.cpp
        TPZPorousElasticResponseTranslator.cpp
        TPZSandlerDimaggioThermoForceATranslator.h
        TPZSandlerDimaggioTranslator.cpp
        TPZSandlerExtendedTranslator.cpp
        TPZYCCamClayPVTranslator.cpp
        TPZYCDruckerPragerPVTranslator.cpp
        TPZYCMohrCoulombPVTranslator.cpp
        TPZYCSandlerDimaggioL2Translator.cpp
        TPZYCSandlerDimaggioLTranslator.cpp
        TPZYCSandlerDimaggioTranslator.cpp
        )
endif()



install(FILES ${headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Material/REAL/Plasticity)

target_sources(pz PRIVATE ${headers} ${sources})