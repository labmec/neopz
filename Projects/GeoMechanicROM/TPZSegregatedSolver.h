//
//  TPZSegregatedSolver.h
//  PZ
//
//  Created by omar on 3/8/17.
//
//

#ifndef TPZSegregatedSolver_h
#define TPZSegregatedSolver_h

#include <stdio.h>
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSimulationData.h"
#include "TPZTransferFunctions.h"
#include "TPZFLuxPressureAnalysis.h"
#include "TPZElasticAnalysis.h"

class TPZSegregatedSolver : public TPZAnalysis {
    
private:
    
    /** @brief define the simulation data */
    TPZSimulationData * fSimulationData;
    
    /** @brief define the elliptic component */
    TPZElasticAnalysis * felliptic;
    
    /** @brief define the parabolic component */
    TPZFLuxPressureAnalysis * fparabolic;
    
    /** @brief define the tranfer object data */
    TPZTransferFunctions * ftransfer;
    
    /** @brief Residue error */
    STATE ferror;
    
    /** @brief Correction variation */
    STATE fdx_norm;
    
    /** @brief number of newton corrections */
    int fk_iterations;
    
public:
    
    /** @brief default constructor  */
    TPZSegregatedSolver();
    
    /** @brief default desconstructor  */
    ~TPZSegregatedSolver();
    
    /** @brief Copy constructor $ */
    TPZSegregatedSolver(const TPZSegregatedSolver &copy);
    
    /** @brief Copy assignemnt operator $ */
    TPZSegregatedSolver &operator=(const TPZSegregatedSolver &other);
    
    
    /** @brief Set the simulation data */
    void SetSimulationData(TPZSimulationData * SimulationData)
    {
        fSimulationData = SimulationData;
    }
    
    /** @brief Get the space generator */
    TPZSimulationData * SimulationData()
    {
        return fSimulationData;
    }
    
    /** @brief Set elliptic */
    void Set_elliptic(TPZElasticAnalysis * elliptic)
    {
        felliptic = elliptic;
    }
    
    /** @brief Get elliptic */
    TPZElasticAnalysis * elliptic()
    {
        return felliptic;
    }
    
    /** @brief Set parabolic */
    void Set_parabolic(TPZFLuxPressureAnalysis * parabolic)
    {
        fparabolic = parabolic;
    }
    
    /** @brief Get parabolic */
    TPZFLuxPressureAnalysis * parabolic()
    {
        return fparabolic;
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
    
    
    
    /** @brief Get k iterations */
    int k_ietrarions(){
        return fk_iterations;
    }
    
    /** @brief Get k iterations */
    void Set_k_ietrarions(int k){
        fk_iterations = k;
    }
    
    /** @brief execute the evolutionary problem */
    void Run_Evolution(std::string elliptic, std::string parabolic);
    
    /** @brief Execute a euler method step */
    void ExcecuteOneStep();
    
    /** @brief Execute a segregated iteration  */
    void SegregatedIteration();
    
    /** @brief PostProcess results */
    void PostProcessStep(std::string elliptic, std::string parabolic);
    
    /** @brief update last state solution */
    void UpdateState();
    
    /** @brief update current state solution */
    void Update_at_n_State();
    
    /** @brief update global state for the new euler step */
    void UpdateGlobalSolution();
    
    /** @brief keep global last state for restart a euler step */
    void KeepGlobalSolution();
    
};



#endif /* TPZSegregatedSolver_h */
