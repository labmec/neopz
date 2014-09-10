//
//  TPZDarcyAnalysis.h
//  PZ
//
//  Created by omar on 9/8/14.
//
//

#ifndef __PZ__TPZDarcyAnalysis__
#define __PZ__TPZDarcyAnalysis__

#include <iostream>
#include "tpzautopointer.h"
#include "TPZFracData.h"
#include "pzcmesh.h"
#include "pzanalysis.h"

/**
 * @author Nathan Shauer and Omar Duran
 * @since 19/08/2014
 * @brief Creates and analyse the data for hydraulic fracturing
 */
class TPZDarcyAnalysis {
    
public:
    
    /** @brief Default Constructor */
    TPZDarcyAnalysis(TPZAutoPointer<TPZFracData> Data);
    
    /** @brief Destructor */
    ~TPZDarcyAnalysis();
    
    /** @brief Runs simulation */
    void Run();
    
    /** @brief Creates geometric mesh */
    TPZGeoMesh * CreateGMesh(const int nel);

    /** @brief Apply a 2D rotation around z axis */
    void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle);

    /** @brief Apply nh times a uniform refinement */
    void UniformRefinement(TPZGeoMesh *gMesh, int nh);
    
    /** @brief Creates fmeshvec[0] for flux with H1 space */
    TPZCompMesh * CreateCMeshFluxHdiv();
    
    /** @brief Creates fmeshvec[1] for pressure with L2 space */
    TPZCompMesh * CreateCMeshPressureL2();
    
    /** @brief Creates Multiphysic mesh for mixed simulation of darcy flow */
    TPZCompMesh * CreateCMeshMixed();
    
    /** @brief Creates interfaces for mixed simulation of darcy flow */
    void CreateInterfaces(TPZCompMesh *cmesh);
    
    /** @brief Assemble last step */
    void AssembleLastStep(TPZAnalysis *an);
    
    /** @brief Newton Method */
    void IterativeProcess(TPZAnalysis *an, std::ostream &out);
    
    /** @brief Solve time steps */
    void SolveSistTransient(TPZAnalysis *an);
    
private:
    
    /** @brief Data of the simulation */
    TPZAutoPointer<TPZFracData> fData;
    
    /** @brief Geometric mesh */
    TPZGeoMesh * fgmesh;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZVec<TPZCompMesh *> fmeshvec;
    
    /** @brief Multphysics cmesh for mixed analysis */
    TPZCompMesh * fcmeshMixed;
    
    /** @brief Mass Residual of step n */
    TPZFMatrix<> fLastStepRhs;
    
};

#endif /* defined(__PZ__TPZDarcyAnalysis__) */