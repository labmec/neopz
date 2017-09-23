//
//  TPZGeomechanicAnalysis.hpp
//  PZ
//
//  Created by Omar on 8/28/16.
//
//

#ifndef TPZGeomechanicAnalysis_hpp
#define TPZGeomechanicAnalysis_hpp

#include <stdio.h>

#include <stdio.h>
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSimulationData.h"

class TPZGeomechanicAnalysis : public TPZAnalysis {
    
private:
    
    /** @brief define the simulation data */
    TPZSimulationData * fSimulationData;
    
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
    
    /** @brief Strain-Stress solution data */
    TPZStack< std::pair<REAL,REAL> > fstrain_stress_duplets;
    
    /** @brief Strain-Porosity solution data */
    TPZStack< std::pair<REAL,REAL> > fstrain_porosity_duplets;
    
    /** @brief Strain-Permeability solution data */
    TPZStack< std::pair<REAL,REAL> > fstrain_permeability_duplets;
    
    /** @brief Strain-Pressure solution data */
    TPZStack< std::pair<REAL,REAL> > fstrain_pressure_duplets;
    
    /** @brief Residue error */
    STATE ferror;
    
    /** @brief Correction variation */
    STATE fdx_norm;
    
    /** @brief number of newton corrections */
    int fk_iterations;
    
public:
    
    /** @brief default constructor  */
    TPZGeomechanicAnalysis();
    
    /** @brief default desconstructor  */
    ~TPZGeomechanicAnalysis();
    
    /** @brief Copy constructor $ */
    TPZGeomechanicAnalysis(const TPZGeomechanicAnalysis &copy);
    
    /** @brief Copy assignemnt operator $ */
    TPZGeomechanicAnalysis &operator=(const TPZGeomechanicAnalysis &other);
    
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
    
    /** @brief PostProcess results for nonlinear case */
    void PostNonlinearProcessStep(std::string plotfile);
    
    /** @brief PostProcess results */
    void PostProcessStep(std::string plotfile);
    
    /** @brief update last state solution */
    void UpdateState();
    
    /** @brief update current state solution */
    void Update_at_n_State();

    /** @brief execute the evolutionary problem */
    void Run_Evolution(TPZVec<REAL> &x, std::string plotfile);
    
    /** @brief Compute the strain and the stress at x euclidean point for each time */
    void AppendStrain_Stress(TPZVec<REAL> & x);
    
    /** @brief Compute the strain and the Porosity at x euclidean point for each time */
    void AppendStrain_Pososity(TPZVec<REAL> & x);
    
    /** @brief Compute the strain and the Permeability at x euclidean point for each time */
    void AppendStrain_Permeability(TPZVec<REAL> & x);
    
    /** @brief Compute the strain and the Pressure at x euclidean point for each time */
    void AppendStrain_Pressure(TPZVec<REAL> & x);
    
    /** @brief Compute the strain and the stress at x euclidean point for each time */
    void PlotStrainStress(std::string file_name);
    
    /** @brief Compute the strain and the Porosity at x euclidean point for each time */
    void PlotStrainPorosity(std::string file_name);
    
    /** @brief Compute the strain and the Permeability at x euclidean point for each time */
    void PlotStrainPermeability(std::string file_name);
    
    /** @brief Compute the strain and the Pressure at x euclidean point for each time */
    void PlotStrainPressure(std::string file_name);
    
};


#endif /* TPZGeomechanicAnalysis_hpp */
