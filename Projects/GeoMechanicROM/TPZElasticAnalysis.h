//
//  TPZElasticAnalysis.h
//  PZ
//
//  Created by Omar on 3/5/17.
//
//

#ifndef TPZElasticAnalysis_h
#define TPZElasticAnalysis_h

#include <stdio.h>
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSimulationData.h"
#include "TPZTransferFunctions.h"


class TPZElasticAnalysis : public TPZAnalysis {
    
private:
    
    /** @brief define the simulation data */
    TPZSimulationData * fSimulationData;
    
    /** @brief define the tranfer object data */
    TPZTransferFunctions * ftransfer;
    
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
    TPZElasticAnalysis();
    
    /** @brief default desconstructor  */
    ~TPZElasticAnalysis();
    
    /** @brief Copy constructor $ */
    TPZElasticAnalysis(const TPZElasticAnalysis &copy);
    
    /** @brief Copy assignemnt operator $ */
    TPZElasticAnalysis &operator=(const TPZElasticAnalysis &other);
    
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
    
    /** @brief Set the simulation data */
    void SetSimulationData(TPZSimulationData * SimulationData)
    {
        fSimulationData = SimulationData;
        fmeshvec.Resize(2);
    }
    
    /** @brief Get the space generator */
    TPZSimulationData * SimulationData()
    {
        return fSimulationData;
    }
    
    /** @brief Set the tranfer object */
    void SetTransfer_object(TPZTransferFunctions * transfer)
    {
        ftransfer = transfer;
    }
    
    /** @brief Get the tranfer object */
    TPZTransferFunctions * Transfer_object()
    {
        return ftransfer;
    }
    
    /** @brief Set vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure */
    void SetMeshvec(TPZVec<TPZCompMesh * > &Meshvec)
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
    
    /** @brief Execute a quasi newton iteration  */
    void QuasiNewtonIteration();
    
    /** @brief Execute a newton iteration  */
    void NewtonIteration();
    
    /** @brief PostProcess results */
    void PostProcessStep(std::string plotfile);
    
    /** @brief update last state solution */
    void UpdateState();
    
    /** @brief update current state solution */
    void Update_at_n_State();

};


#endif /* TPZElasticAnalysis_h */
