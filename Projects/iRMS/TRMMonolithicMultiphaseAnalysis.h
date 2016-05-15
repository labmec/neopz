//
//  TRMMonolithicMultiphaseAnalysis.h
//  PZ
//
//  Created by Omar on 5/11/16.
//
//

#ifndef TRMMonolithicMultiphaseAnalysis_h
#define TRMMonolithicMultiphaseAnalysis_h

#include <stdio.h>
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TRMSpaceOdissey.h"
#include "TRMSimulationData.h"


class TRMMonolithicMultiphaseAnalysis : public TPZAnalysis {
    
private:
    
    /** @brief define the simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;

    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2, fmeshvec[2] = SaturationL2 */
    TPZManVector<TPZCompMesh * , 4> fmeshvec;

    /** @brief Part of residue at n state  */
    TPZFMatrix<STATE> fResidue_n;
    
    /** @brief Part of residue at past state  */
    TPZFMatrix<STATE> fResidue;
    
    /** @brief Solution ate n state */
    TPZFMatrix<STATE> fSolution_n;
    
    /** @brief Solution at past state */
    TPZFMatrix<STATE> fSolution;
    
    /** @brief Residue error */
    STATE ferror;
    
    /** @brief Correction variation */
    STATE fdx_norm;
    
public:
    
    /** @brief default constructor  */
    TRMMonolithicMultiphaseAnalysis();

    /** @brief default desconstructor  */    
    ~TRMMonolithicMultiphaseAnalysis();
    
    
    /**
     * @defgroup Access Methods
     * @brief    Implements Access methods:
     * @{
     */
    
    /** @brief Set the simulation data */
    void SetSimulationData(TPZAutoPointer<TRMSimulationData> SimulationData)
    {
        fSimulationData = SimulationData;
    }
    /** @brief Get the space generator */
    TPZAutoPointer<TRMSimulationData> SimulationData()
    {
        return fSimulationData;
    }
    
    // @}
    
    /**
     * @defgroup Time foward methods
     * @{
     */
    
    /** @brief Execute a euler method step */
    void ExcecuteOneStep();

    /** @brief Execute a newton iteration  */
    void NewtonIteration();
    
    /** @brief PostProcess results */
    void PostProcessStep();
    
    // @}
    

    
    
};

#endif /* MonolithicMultiphaseAnalysis_h */
