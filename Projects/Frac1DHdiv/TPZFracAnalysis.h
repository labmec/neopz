#ifndef TPZFracAnalysis_H
#define TPZFracAnalysis_H

#include "tpzautopointer.h"
#include "TPZFracData.h"
#include "pzcmesh.h"
#include "pzanalysis.h"

class TPZMatfrac1dhdiv;

/**
 * @author Omar Duran and Nathan Shauer
 * @since 19/08/2014
 * @brief Creates and analyse the data for hydraulic fracturing
 */
class TPZFracAnalysis {
  
public:
  
  /** @brief Default Constructor */
  TPZFracAnalysis(TPZAutoPointer<TPZFracData> Data);
  
  /** @brief Destructor */
  ~TPZFracAnalysis();
  
  /** @brief Runs case where fracture starts from zero */
  void Run();
  
  /** @brief Runs test where fracture has determined length */
  void RunTest();
  
  /** @brief Creates geometric mesh */
  TPZGeoMesh * CreateGMesh();
  
  /** @brief Creates fmeshvec[0] for flux with H1 space */
  TPZCompMesh * CreateCMeshFluxH1();
  
  /** @brief Creates fmeshvec[1] for pressure with L2 space */
  TPZCompMesh * CreateCMeshPressureL2();
  
  /** @brief Creates Multiphysic mesh for mixed simulation of fracture */
  TPZCompMesh * CreateCMeshMixed(TPZFMatrix<REAL> vlMatrix);
  
  /** @brief Assemble last step */
  void AssembleLastStep(TPZAnalysis *an);
  
  /** @brief Newton Method */
  void IterativeProcess(TPZAnalysis *an, std::ostream &out, int nit = 50);
  
  /** @brief Solve time steps */
  bool SolveSistTransient(TPZAnalysis *an);
  
  /** @brief Updates Leak Off integration points values */
  void AcceptSolution(TPZAnalysis *an);
  
  /** @brief PostProcess mesh in VTK */
  void PostProcessVTK(TPZAnalysis *an);
  
  /** @brief Calculates Q of the tip of the fracture */
  REAL Qtip();
  
  /** @brief Set the pressure on last element as equal of the second last element */
  void SetPressureOnLastElement(TPZAnalysis *an);
  
  /** @brief Return the quantity of fracture elements */
  int HowManyFracElement();
  
  /** @brief Initialize x0 for newton iterations */
  void ComputeFirstSolForOneELement(TPZAnalysis * an);
  
  /** @brief Finds the initial time step to run simulation and returns the vl to train the integration points */
  REAL RunUntilOpen();
  
  /** @brief Return the Q (flow) that a fresh new element could absorb in one time step */
  REAL QOfAFreshNewElement();
  
  /** @brief Creates first GeoEl with bc */
  TPZGeoEl* CreateFirstGeoElWithBC();
  
  /** @brief Verifies if has to propagate, ie, the qtip is bigger than leak off of the next element */
  bool VerifyIfPropagate(REAL qtip);
  
  /** @brief Find the pressure BC geo element */
  TPZGeoEl * FindPressureBCElement();
  
  /** @brief Return the flow criteria used to decide if propagates fracture */
  REAL PropagationFlowCriteria(REAL qFreshNewEl, REAL ql);
  
private:
  
  /** @brief bool which indicates if the end of time is reached */
  bool fmustStop;
  
  /** @brief Data of the simulation */
  TPZAutoPointer<TPZFracData> fData;
  
  /** @brief Geometric mesh */
  TPZGeoMesh * fgmesh;
  
  /** @brief Vector of compmesh pointers. fmeshvec[0] = flowH1, fmeshvec[1] = PressureL2 */
  TPZVec<TPZCompMesh *> fmeshvec;
  
  /** @brief Multphysics cmesh for mixed analysis */
  TPZCompMesh * fcmeshMixed;
  
  /** @brief Mass Residual of step n */
  TPZFMatrix<STATE> fLastStepRhs;
  
  /** @brief Pointer to material of fracturing phenomena  */
  TPZMatfrac1dhdiv *fMatFrac;
  
};

#endif