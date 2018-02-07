#ifndef COMMON3DHPP
#define COMMON3DHPP

#include "pzanalysis.h"
#include "pzcmesh.h"
#include "TPZAnalyticSolution.h"

#ifdef _AUTODIFF
extern TLaplaceExample1 ExactLaplace;

extern TElasticity3DAnalytic ExactElast;
#endif

//    This Solve Different analysis
void SolveSist(TPZAnalysis *an, TPZCompMesh *fCmesh, int numthreads);

/// insert material objects in the computational mesh
void InsertMaterialObjects3D(TPZCompMesh *cmesh, bool scalarproblem);

/// Build a square mesh with boundary conditions
TPZCompMesh *SetupSquareMesh3D(int nelx, int nrefskeleton, int porder, bool elasticityproblem);

enum MMATID {Ebc1 = -1,Enomat, Emat1, Emat2, Emat3, Emat4, Ebc2, Ebc3, Ebc4, Ebc5, Ebcpoint1, Ebcpoint2, Ebcpoint3, Ewrap, ESkeleton, EInterfaceMat1, EInterfaceMat2, EGroup};

#ifdef _AUTODIFF
/// Function defining the exact elasticity solution
inline void Elasticity_exact(const TPZVec<REAL> &xv, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv)
{
    ExactElast.Solution(xv,val,deriv);
}
#endif

/// Read a UNSWSBFem file
TPZGeoMesh *ReadUNSWSBGeoFile(const std::string &filename, TPZVec<int64_t> &elpartition, TPZVec<int64_t> &scalingcenterindices);

#endif
