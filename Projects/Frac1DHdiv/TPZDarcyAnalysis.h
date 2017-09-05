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

    /** @brief Creates a L2 projection for intial pressure */
    TPZCompMesh * L2ProjectionP(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);
    
    /** @brief Creates Multiphysic mesh for mixed simulation of darcy flow */
    TPZCompMesh * CreateCMeshMixed(TPZFMatrix<REAL> Vl);
    
    /** @brief Creates interfaces for mixed simulation of darcy flow */
    void CreateInterfaces();
  
  
    
    /** @brief Assemble last step */
    void AssembleLastStep(TPZAnalysis *an);
    
    /** @brief Newton Method -> Solves  X(u) [DeltaU] = - [R(u)]; */
    void IterativeProcess(TPZAnalysis *an, std::ostream &out, int numiter = 50);
    
    /** @brief Solve LS */
    void SolveSyst(TPZAnalysis &an, TPZCompMesh *Cmesh);
    
    /** @brief Solve time steps */
    bool SolveSistTransient(TPZAnalysis *an, bool initial);

    /** @brief Solve time steps with Fracture elements */
    bool SolveSistTransientWithFracture(TPZAnalysis *an);
    
    /** @brief PostProcess mesh in VTK */
    void PostProcessVTK(TPZAnalysis *an);
    
    /** @brief Print geometric mesh in VTK */
    void PrintGeometricMesh(TPZGeoMesh * gmesh);

    /** @brief Print geometric mesh in VTK */
    void PrintComputatiolaMeshInfo(TPZCompMesh * cmesh);

    /** @brief adjust mtid's for propagation */
    void InsertFracGeoMesh();
  
    /** @brief adjust cmesh */
    void InsertFracCompMesh();

    /** @brief sets the memory of the new integration points */
    void SetIntPointMemory();
  
    /** @brief Creates a quadrilateral conform gmesh with max_edge defined by param */
    void CreateGeoMeshQuad(int npropag, int nrefy, REAL ly);
  
    /** @brief adjust mtid's for propagation */
    void DarcyGmesh(TPZGeoMesh * gmesh);
    
    /** @brief Creat Computational Multiphysics Mesh */
    void CreateMultiphysicsMesh(TPZFMatrix<REAL> Vl);

    /** @brief Updates Leak Off integration points values */
    void AcceptSolution(TPZAnalysis *an);

    /** @brief Set the pressure on last element as equal of the second last element */
    void SetPressureOnLastElement(TPZAnalysis *an);

    /** @brief Set the pressure on a new element in order to get zero opening for assemble last step */
    void SetPressureOnNewElement(TPZAnalysis *an);
    
    /** @brief Initialize x0 for newton iterations */    
    void ComputeFirstSolForOneELement(TPZAnalysis * an);
    
    /** @brief Return the quantity of fracture elements */
    int HowManyFracElement();
    
    /** @brief Calculates Q of the tip of the fracture */
    REAL Qtip();
  
    /** @brief Finds the initial time step to run simulation and returns the vl to train the integration points */
    REAL RunUntilOpen();
    
    /** @brief Return the Q (flow) that a fresh new element could absorb in one time step */
    REAL QOfAFreshNewElement();
    
    /** @brief Creates first GeoEl with bc */
    TPZGeoEl* CreateFirstGeoElWithBC();
    
    /** @brief Verifies if has to propagate, ie, the qtip is bigger than leak off of the next element */
    bool VerifyIfPropagate(REAL qtip);
    
    /** @brief Return the flow criteria used to decide if propagates fracture */
    REAL PropagationFlowCriteria(REAL qFreshNewEl, REAL ql);
  
    /** @brief switch hdivboud element to discontinuous galerkin element */
    void SwitchTipElement(TPZCompEl * cel, TPZCompEl *celpoint, TPZGeoEl *gelpoint);
  
    /** @brief Sets the fconnectvec of interface multiphysics elements */
    void SetInterfaceConnects();
  
    /** @brief Switch hdivboud element to auxiliar qtip inlet */
    void SwitchBCInFrontOfFrac(TPZGeoEl * gel);

  
private:
    
    /** @brief bool which indicates if the end of time is reached */
    bool fmustStop;
    
    /** @brief Data of the simulation */
    TPZAutoPointer<TPZFracData> fData;
    
    /** @brief Geometric mesh */
    TPZGeoMesh * fgmesh;
    
    /** @brief Vector of compmesh pointers. fmeshvec[0] = flowHdiv, fmeshvec[1] = PressureL2 */
    TPZManVector<TPZCompMesh * , 2> fmeshvec;
    
    /** @brief Multphysics cmesh for mixed analysis */
    TPZCompMesh * fcmeshMixed;
    
    /** @brief Mass Residual of step n */
    TPZFMatrix<STATE> fLastStepRhs;
    
};

#endif /* defined(__PZ__TPZDarcyAnalysis__) */
