//
//  TRMGeomechanicAnalysis.h
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#ifndef __PZ__TRMGeomechanicAnalysis__
#define __PZ__TRMGeomechanicAnalysis__

#include <stdio.h>
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TRMSimulationData.h"
#include "TRMBuildTransfers.h"

class TRMGeomechanicAnalysis : public TPZAnalysis{
    
private:
    
    /** @brief define the simulation data */
    TRMSimulationData * fSimulationData;
    
    /** @brief define the transfer matrices */
    TRMBuildTransfers * fTransfer;
    
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
    
    /** @brief number of newton corrections */
    int fk_iterations;
    
public:
    
    /** @brief default constructor  */
    TRMGeomechanicAnalysis();
    
    /** @brief default desconstructor  */
    ~TRMGeomechanicAnalysis();
    
    /** @brief Copy constructor $ */
    TRMGeomechanicAnalysis(const TRMGeomechanicAnalysis &copy);
    
    /** @brief Copy assignemnt operator $ */
    TRMGeomechanicAnalysis &operator=(const TRMGeomechanicAnalysis &other);
    
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
    
    
};

#endif /* defined(__PZ__TRMGeomechanicsAnalysis__) */
