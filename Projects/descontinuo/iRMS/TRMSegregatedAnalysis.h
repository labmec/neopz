//
//  TRMSegregatedAnalysis.hpp
//  PZ
//
//  Created by Omar on 7/18/16.
//
//

#ifndef TRMSegregatedAnalysis_h
#define TRMSegregatedAnalysis_h

#include <stdio.h>

#include <stdio.h>
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TRMSimulationData.h"
#include "TRMBuildTransfers.h"
#include "TRMFluxPressureAnalysis.h"
#include "TRMTransportAnalysis.h"

class TRMSegregatedAnalysis : public TPZAnalysis {
    
private:
    
    /** @brief define the simulation data */
    TRMSimulationData * fSimulationData;
    
    /** @brief define the transfer matrices */
    TRMBuildTransfers * fTransfer;

    /** @brief define the parabolic system */
    TRMFluxPressureAnalysis * fParabolic_mhm;
    
    /** @brief define the parabolic system */
    TRMFluxPressureAnalysis * fParabolic;
    
    /** @brief define the hyperbolic system */
    TRMTransportAnalysis * fHyperbolic;
    
    /** @brief Residue error for flux - pressure */
    REAL ferror_flux_pressure;
    
    /** @brief Residue error for saturations */
    REAL ferror_saturation;
    
    /** @brief Correction variation for flux - pressure */
    REAL fdx_norm_flux_pressure;
    
    /** @brief Correction variation for saturations */
    REAL fdx_norm_saturation;

    
public:
    
    /** @brief default constructor  */
    TRMSegregatedAnalysis();
    
    /** @brief default desconstructor  */
    ~TRMSegregatedAnalysis();
    
    /** @brief Copy constructor $ */
    TRMSegregatedAnalysis(const TRMSegregatedAnalysis &copy);
    
    /** @brief Copy assignemnt operator $ */
    TRMSegregatedAnalysis &operator=(const TRMSegregatedAnalysis &other);
    
    /**
     * @defgroup Access Methods
     * @brief    Implements Access methods:
     * @{
     */
    
    
    /** @brief Set the simulation data */
    void SetSimulationData(TRMSimulationData * SimulationData)
    {
        fSimulationData = SimulationData;
    }
    
    /** @brief Get the space generator */
    TRMSimulationData * SimulationData()
    {
        return fSimulationData;
    }
    
    /** @brief Set the transfer object */
    void SetTransfer(TRMBuildTransfers *Transfer)
    {
        fTransfer = Transfer;
    }
    
    /** @brief Get the transfer object */
    TRMBuildTransfers * Transfer()
    {
        return fTransfer;
    }
    
    /** @brief Set the parabolic mhm part */
    void SetParabolicMHM(TRMFluxPressureAnalysis * parabolic_mhm)
    {
        fParabolic_mhm = parabolic_mhm;
    }
    
    /** @brief Get the parabolic mhm part  */
    TRMFluxPressureAnalysis * ParabolicMHM()
    {
        return fParabolic_mhm;
    }

    /** @brief Set the parabolic part */
    void SetParabolic(TRMFluxPressureAnalysis * parabolic)
    {
        fParabolic = parabolic;
    }
    
    /** @brief Get the parabolic part  */
    TRMFluxPressureAnalysis * Parabolic()
    {
        return fParabolic;
    }
    
    /** @brief Set the hyperbolic part */
    void SetHyperbolic(TRMTransportAnalysis * hyperbolic)
    {
        fHyperbolic = hyperbolic;
    }
    
    /** @brief Get the hyperbolic part  */
    TRMTransportAnalysis * Hyperbolic()
    {
        return fHyperbolic;
    }
    
    /** @brief Resize and fill residue and solution vectors */
    void AdjustVectors();
    
    // @}
    
    /**
     * @defgroup Time foward methods
     * @{
     */
    
    /** @brief Execute a euler method step */
    void ExcecuteOneStep();
    
    /** @brief Execute a newton iteration  */
    void NewtonIteration();
    
    /** @brief Execute a segregated iteration  */
    void SegregatedIteration();
    
    /** @brief PostProcess results */
    void PostProcessStep(bool draw_mixed_mapQ);

    /** @brief Update memory using the Transfer object at REAL n */
    void UpdateFluxes_at_n();
    
    /** @brief Update memory using the Transfer object at REAL n */
    void UpdateMemory_at_n();
    
    /** @brief Update memory using the Transfer object */
    void UpdateMemory();
    
    /** @brief update global state for the new euler step */
    void UpdateGlobalSolution();
    
    /** @brief keep global last state for restart a euler step */
    void KeepGlobalSolution();
    
    // @}
    
    
    
    
};


#endif /* TRMSegregatedAnalysis_h */
