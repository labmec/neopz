//
//  TRMTransportAnalysis.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMTransportAnalysis__
#define __PZ__TRMTransportAnalysis__

#include <stdio.h>
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TRMSimulationData.h"
#include "TRMBuildTransfers.h"

class TRMTransportAnalysis : public TPZAnalysis{
        
private:
    
    /** @brief define the simulation data */
    TPZAutoPointer<TRMSimulationData> fSimulationData;
    
    /** @brief define the transfer matrices */
    TPZAutoPointer<TRMBuildTransfers> fTransfer;
    
    /** @brief Vector of compmesh pointers. fmeshvec = { alpha, beta, gamma , etc...} phases */
    TPZVec<TPZCompMesh *> fmeshvec;
    
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
    TRMTransportAnalysis();
    
    /** @brief default desconstructor  */
    ~TRMTransportAnalysis();
    
    /** @brief Copy constructor $ */
    TRMTransportAnalysis(const TRMTransportAnalysis &copy);
    
    /** @brief Assignemnt operator $ */
    TRMTransportAnalysis &operator=(const TRMTransportAnalysis &other);
    
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
    
    /** @brief Set vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure */
    void SetMeshvec(TPZVec<TPZCompMesh *> &Meshvec)
    {
        fmeshvec = Meshvec;
    }
    /** @brief Get Vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure */
    TPZVec<TPZCompMesh *> & Meshvec()
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

#endif /* defined(__PZ__TRMTransportAnalysis__) */