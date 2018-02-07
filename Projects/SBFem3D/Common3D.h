#ifndef COMMON3DHPP
#define COMMON3DHPP

#include "pzanalysis.h"
#include "pzcmesh.h"
#include "TPZAnalyticSolution.h"

//    This Solve Different analysis
void SolveSist(TPZAnalysis *an, TPZCompMesh *fCmesh);

/// insert material objects in the computational mesh
void InsertMaterialObjects3D(TPZCompMesh *cmesh, int problemtype);

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupSquareMesh3D(int nelx, int nrefskeleton, int porder, bool elasticityproblem);

enum MMATID {Ebc1 = -1,Enomat, Emat1, Emat2, Emat3, Emat4, Ebc2, Ebc3, Ebc4, Ewrap, ESkeleton, EInterfaceMat1, EInterfaceMat2, EGroup};

/// Function defining the Harmonic solution at the left of the domain
void HarmonicNeumannLeft(const TPZVec<REAL> &x, TPZVec<STATE> &val);

/// Function defining the Harmonic solution at the right of the domain
void HarmonicNeumannRight(const TPZVec<REAL> &x, TPZVec<STATE> &val);

/// Function defining the exact harmonic solution
void Harmonic_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv);

extern TLaplaceExampleSmooth ExactLaplace;

#endif
