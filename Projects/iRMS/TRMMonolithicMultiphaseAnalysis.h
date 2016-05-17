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
    TPZFMatrix<STATE> fR_n;
    
    /** @brief Part of residue at past state  */
    TPZFMatrix<STATE> fR;
    
    /** @brief Solution ate n state */
    TPZFMatrix<STATE> fX_n;
    
    /** @brief Solution at past state */
    TPZFMatrix<STATE> fX;
    
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
    void SetSimulationData(TPZAutoPointer<TRMSimulationData> &SimulationData)
    {
        fSimulationData = SimulationData;
        if (fSimulationData->IsOnePhaseQ()) {
            fmeshvec.Resize(2);
        }
        
        if (fSimulationData->IsTwoPhaseQ()) {
            fmeshvec.Resize(3);
        }
        
        if (fSimulationData->IsThreePhaseQ()) {
            fmeshvec.Resize(4);
        }
        
    }
    
    /** @brief Get the space generator */
    TPZAutoPointer<TRMSimulationData> SimulationData()
    {
        return fSimulationData;
    }
    
    /** @brief Set vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure, fmeshvec[2] = Water, fmeshvec[3] = Oil */
    void SetMeshvec(TPZManVector<TPZCompMesh * , 4> &Meshvec)
    {
        fmeshvec = Meshvec;
    }
    /** @brief Get Vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure, fmeshvec[2] = Water, fmeshvec[3] = Oil */
    TPZManVector<TPZCompMesh * , 4> & Meshvec()
    {
        return fmeshvec;
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
    
    /** @brief PostProcess results */
    void PostProcessStep();
    
    // @}
    

    
    
};

#endif /* MonolithicMultiphaseAnalysis_h */
