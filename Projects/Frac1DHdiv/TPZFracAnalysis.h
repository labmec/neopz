#ifndef TPZFracAnalysis_H
#define TPZFracAnalysis_H

#include "tpzautopointer.h"
#include "TPZFracData.h"
#include "pzcmesh.h"
#include "pzanalysis.h"

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

  /** @brief Runs simulation */
  void Run();
  
  /** @brief Creates geometric mesh */
  TPZGeoMesh * CreateGMesh(const int nel);
  
  /** @brief Creates fmeshvec[0] for flux with H1 space */
  TPZCompMesh * CreateCMeshFluxH1();
  
  /** @brief Creates fmeshvec[1] for pressure with L2 space */
  TPZCompMesh * CreateCMeshPressureL2();
  
  /** @brief Creates Multiphysic mesh for mixed simulation of fracture */
  TPZCompMesh * CreateCMeshMixed();
  
  /** @brief Newton Method */
  void IterativeProcess(TPZAnalysis *an, std::ostream &out);
  
private:
  
  /** @brief Data of the simulation */
  TPZAutoPointer<TPZFracData> fData;
  
  /** @brief Geometric mesh */
  TPZGeoMesh * fgmesh;
  
  /** @brief Vector of compmesh pointers. fmeshvec[0] = flowH1, fmeshvec[1] = PressureL2 */
  TPZVec<TPZCompMesh *> fmeshvec;
  
  /** @brief Multphysics cmesh for mixed analysis */
  TPZCompMesh * fcmeshMixed;
  
};

#endif