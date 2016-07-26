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
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    /** @brief define the transfer matrices */
    TPZAutoPointer<TRMBuildTransfers> fTransfer;
    
    /** @brief define the parabolic system */
    TPZAutoPointer<TRMFluxPressureAnalysis> fParabolic;
    
    /** @brief define the hyperbolic system */
    TPZAutoPointer<TRMTransportAnalysis> fHyperbolic;
    
    /** @brief Residue error for flux - pressure */
    STATE ferror_flux_pressure;
    
    /** @brief Residue error for saturations */
    STATE ferror_saturation;
    
    /** @brief Correction variation for flux - pressure */
    STATE fdx_norm_flux_pressure;
    
    /** @brief Correction variation for saturations */
    STATE fdx_norm_saturation;

    
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
    void SetSimulationData(TPZAutoPointer<TRMSimulationData> &SimulationData)
    {
        fSimulationData = SimulationData;
    }
    
    /** @brief Get the space generator */
    TPZAutoPointer<TRMSimulationData> SimulationData()
    {
        return fSimulationData;
    }
    
    /** @brief Set the transfer object */
    void SetTransfer(TPZAutoPointer<TRMBuildTransfers> &Transfer)
    {
        fTransfer = Transfer;
    }
    
    /** @brief Get the transfer object */
    TPZAutoPointer<TRMBuildTransfers> Transfer()
    {
        return fTransfer;
    }

    /** @brief Set the parabolic part */
    void SetParabolic(TPZAutoPointer<TRMFluxPressureAnalysis> &parabolic)
    {
        fParabolic = parabolic;
    }
    
    /** @brief Get the parabolic part  */
    TPZAutoPointer<TRMFluxPressureAnalysis> Parabolic()
    {
        return fParabolic;
    }
    
    /** @brief Set the hyperbolic part */
    void SetHyperbolic(TPZAutoPointer<TRMTransportAnalysis> &hyperbolic)
    {
        fHyperbolic = hyperbolic;
    }
    
    /** @brief Get the hyperbolic part  */
    TPZAutoPointer<TRMTransportAnalysis> Hyperbolic()
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
    void ExcecuteOneStep(bool flag);
    
    /** @brief Execute a newton iteration  */
    void NewtonIteration();
    
    /** @brief Execute a segregated iteration  */
    void SegregatedIteration();
    
    /** @brief PostProcess results */
    void PostProcessStep();

    /** @brief Update memory using the Transfer object at state n */
    void UpdateFluxes_at_n();
    
    /** @brief Update memory using the Transfer object at state n */
    void UpdateMemory_at_n();
    
    /** @brief Update memory using the Transfer object */
    void UpdateMemory();
    
    
    // @}
    
    
    
    
};


#endif /* TRMSegregatedAnalysis_h */
