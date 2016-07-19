//
//  TRMFluxPressureAnalysis.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMFluxPressureAnalysis__
#define __PZ__TRMFluxPressureAnalysis__

#include <stdio.h>
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TRMSimulationData.h"
#include "TRMBuildTransfers.h"

class TRMFluxPressureAnalysis : public TPZAnalysis {
    
private:
    
    /** @brief define the simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    /** @brief define the transfer matrices */
    TPZAutoPointer<TRMBuildTransfers> fTransfer;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 2> fmeshvec;
    
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
    TRMFluxPressureAnalysis();
    
    /** @brief default desconstructor  */
    ~TRMFluxPressureAnalysis();
    
    /** @brief Copy constructor $ */
    TRMFluxPressureAnalysis(const TRMFluxPressureAnalysis &copy);
    
    /** @brief Copy assignemnt operator $ */
    TRMFluxPressureAnalysis &operator=(const TRMFluxPressureAnalysis &other);
    
    /**
     * @defgroup Access Methods
     * @brief    Implements Access methods:
     * @{
     */
    
    /** @brief Set Solution at n state */
    void SetX_n(TPZFMatrix<STATE> &x){
        fX_n = x;
    }
    
    /** @brief Set Solution at n state */
    TPZFMatrix<STATE> & X_n(){
        return fX_n;
    }
    
    /** @brief Set Solution at past state */
    void SetX(TPZFMatrix<STATE> &x){
        fX = x;
    }
    
    /** @brief Set Solution at past state */
    TPZFMatrix<STATE> & X(){
        return fX;
    }
    
    /** @brief Set the simulation data */
    void SetSimulationData(TPZAutoPointer<TRMSimulationData> &SimulationData)
    {
        fSimulationData = SimulationData;
        fmeshvec.Resize(2);
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
    
    /** @brief Set vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure */
    void SetMeshvec(TPZManVector<TPZCompMesh * , 2> &Meshvec)
    {
        fmeshvec = Meshvec;
    }
    /** @brief Get Vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure */
    TPZManVector<TPZCompMesh * , 2> & Meshvec()
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
    
    /** @brief Update memory using the Transfer object at state n */
    void UpdateMemory_at_n();

    /** @brief Update memory using the Transfer object */
    void UpdateMemory();
    
    
    // @}

    
    
    
};

#endif /* defined(__PZ__TRMFluxPressureAnalysis__) */
