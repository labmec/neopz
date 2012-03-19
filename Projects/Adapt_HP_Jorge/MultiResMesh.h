//$Id: MultiResMesh.h,v 1.2 2009-11-04 14:13:24 fortiago Exp $

class TPZGeoMesh;
class TPZCompMesh;
class TPZGeoEl;
#include <set>
#include <pzvec.h>

template <class TVar>
class TPZFMatrix;


/** Principal work */
TPZCompMesh *CreateMeshMultires(TPZGeoMesh * gmesh);
void InitialSolutionMultires(TPZFMatrix<REAL> &InitialSol, TPZCompMesh * cmesh);

double ComputeTimeStep(double CFL, int Level, TPZGeoMesh * gmesh);

TPZGeoMesh * CreateCoarseMesh(int nLevel);
