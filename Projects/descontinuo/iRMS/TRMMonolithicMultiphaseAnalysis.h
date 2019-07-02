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
    TRMSimulationData * fSimulationData;

    /** @brief Vector of compmesh pointers. fmeshvec[0] = Biot, fmeshvec[1] = flowHdiv, fmeshvec[2] = PressureL2, fmeshvec[3] = SaturationL2 */
    TPZManVector<TPZCompMesh * , 5> fmeshvec;

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
    TRMMonolithicMultiphaseAnalysis();

    /** @brief default desconstructor  */    
    ~TRMMonolithicMultiphaseAnalysis();
    
    /** @brief default constructor  */
    TRMMonolithicMultiphaseAnalysis(const TRMMonolithicMultiphaseAnalysis &copy)
    {
        DebugStop();
    }
    
    /** @brief default constructor  */
    TRMMonolithicMultiphaseAnalysis &operator=(const TRMMonolithicMultiphaseAnalysis &copy)
    {
        DebugStop();
        return *this;
    }
    
    /**
     * @defgroup Access Methods
     * @brief    Implements Access methods:
     * @{
     */
    
    /** @brief Set Solution at n state */
    void SetX_n(TPZFMatrix<STATE> x){
        fX_n = x;
    }
    
    /** @brief Set Solution at n state */
    TPZFMatrix<STATE> X_n(){
        return fX_n;
    }

    /** @brief Set Solution at past state */
    void SetX(TPZFMatrix<STATE> x){
        fX = x;
    }
    
    /** @brief Set Solution at past state */
    TPZFMatrix<STATE> X(){
        return fX;
    }
    
    
    /** @brief Set the simulation data */
    void SetSimulationData(TRMSimulationData * SimulationData)
    {
        fSimulationData = SimulationData;
        if (fSimulationData->IsOnePhaseQ()) {
            fmeshvec.Resize(3);
        }
        
        if (fSimulationData->IsTwoPhaseQ()) {
            fmeshvec.Resize(4);
        }
        
        if (fSimulationData->IsThreePhaseQ()) {
            fmeshvec.Resize(5);
        }
        
    }
    
    /** @brief Get the space generator */
    TRMSimulationData * SimulationData()
    {
        return fSimulationData;
    }
    
    /** @brief Set vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure, fmeshvec[2] = Water, fmeshvec[3] = Oil */
    void SetMeshvec(TPZManVector<TPZCompMesh * , 5> &Meshvec)
    {
        fmeshvec = Meshvec;
    }
    /** @brief Get Vector of compmesh pointers. fmeshvec[0] = flux, fmeshvec[1] = Pressure, fmeshvec[2] = Water, fmeshvec[3] = Oil */
    TPZManVector<TPZCompMesh * , 5> & Meshvec()
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
    
    // @}
    

    
    
};

#endif /* MonolithicMultiphaseAnalysis_h */
