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
    TRMSimulationData * fSimulationData;
    
    /** @brief define the transfer matrices */
    TRMBuildTransfers * fTransfer;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 2> fmeshvec;
    
    /** @brief Part of residue at n REAL  */
    TPZFMatrix<REAL> fR_n;
    
    /** @brief Part of residue at past REAL  */
    TPZFMatrix<REAL> fR;
    
    /** @brief Solution ate n REAL */
    TPZFMatrix<REAL> fX_n;
    
    /** @brief Solution at past REAL */
    TPZFMatrix<REAL> fX;
    
    /** @brief Residue error */
    REAL ferror;
    
    /** @brief Correction variation */
    REAL fdx_norm;
    
    /** @brief number of newton corrections */
    int fk_iterations;
    
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
    
    /** @brief Set Solution at n REAL */
    void SetX_n(TPZFMatrix<REAL> &x){
        fX_n = x;
    }
    
    /** @brief Set Solution at n REAL */
    TPZFMatrix<REAL> & X_n(){
        return fX_n;
    }
    
    /** @brief Set Solution at past REAL */
    void SetX(TPZFMatrix<REAL> &x){
        fX = x;
    }
    
    /** @brief Set Solution at past REAL */
    TPZFMatrix<REAL> & X(){
        return fX;
    }
    
    /** @brief Set the simulation data */
    void SetSimulationData(TRMSimulationData * SimulationData)
    {
        fSimulationData = SimulationData;
        fmeshvec.Resize(2);
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
    
    // @}
    
    /**
     * @defgroup Time foward methods
     * @{
     */
    
    /** @brief Execute a euler method step */
    void ExcecuteOneStep();
    
    /** @brief Execute a newton iteration  */
    void NewtonIteration();

    /** @brief Execute a quasi newton iteration  */
    void QuasiNewtonIteration();
    
    /** @brief PostProcess results */
    void PostProcessStep();
    
    /** @brief Update memory using the Transfer object at REAL n */
    void UpdateMemory_at_n();

    /** @brief Update memory using the Transfer object */
    void UpdateMemory();
    
    
    // @}

    
    
    
};

#endif /* defined(__PZ__TRMFluxPressureAnalysis__) */
