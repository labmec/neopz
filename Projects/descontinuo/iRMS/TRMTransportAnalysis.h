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
#include "pzgradientreconstruction.h"

class TRMTransportAnalysis : public TPZAnalysis{
        
private:
    
    /** @brief define the simulation data */
    TRMSimulationData * fSimulationData;
    
    /** @brief define the transfer matrices */
    TRMBuildTransfers * fTransfer;
    
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
    
    /** @brief number of newton corrections */
    int fk_iterations;
    
    /** @brief active equations of alpha saturation */
    TPZManVector<int64_t> factive_sa;
    
    /** @brief no active equations of alpha saturation */
    TPZManVector<int64_t> fno_active_sa;
    
    /** @brief active equations of beta saturation */
    TPZManVector<int64_t> factive_sb;
    
    /** @brief no active equations of beta saturation */
    TPZManVector<int64_t> fno_active_sb;
    
    /** @brief Gradient reconstruction object */
    TPZGradientReconstruction *fgradreconst;
    
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
    void SetTransfer(TRMBuildTransfers * Transfer)
    {
        fTransfer = Transfer;
    }
    
    /** @brief Get the transfer object */
    TRMBuildTransfers * Transfer()
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
    
    /** @brief Get current error */
    REAL error_norm(){
        return ferror;
    }
    
    /** @brief Set dx error */
    void Set_error_norm(REAL error){
        ferror = error;
    }
    
    /** @brief Get dx error */
    REAL dx_norm(){
        return fdx_norm;
    }
    
    /** @brief Set current error */
    void Set_dx_norm(REAL dx_norm){
        fdx_norm = dx_norm;
    }
    
    /** @brief Get k iterations */
    int k_ietrarions(){
        return fk_iterations;
    }
    
    /** @brief Get k iterations */
    void Set_k_ietrarions(int k){
        fk_iterations = k;
    }
    
    
    /** @brief Execute a euler method step */
    void ExcecuteOneStep();
    
    /** @brief Execute a newton iteration  */
    void NewtonIteration();
    
    /** @brief Execute a quasi newton iteration  */
    void QuasiNewtonIteration();
    
    /** @brief PostProcess results */
    void PostProcessStep();
    
    /** @brief Update memory using the Transfer object at state n */
    void UpdateMemory_at_n();
    
    /** @brief Update memory using the Transfer object */
    void UpdateMemory();

    /** @brief Filter Saturations on the non linear system */
    void FilterSaturationGradients();
    
    /** @brief clean up gradiens on solution vector */
    void CleanUpGradients();
    
    /** @brief Saturation reconstruction */
    void SaturationReconstruction();
    
    /** @brief Filtering equations for saturation reconstruction */
    void FilterEquations();
        
    
    
    
};

#endif /* defined(__PZ__TRMTransportAnalysis__) */