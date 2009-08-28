//$Id: MultiResMesh.h,v 1.1 2009-08-28 22:59:35 fortiago Exp $

class TPZGeoMesh;
class TPZCompMesh;
class TPZGeoEl;
#include <set>
#include <pzvec.h>
class TPZFMatrix;


/** Principal work */
TPZCompMesh *CreateMeshMultires(TPZGeoMesh * gmesh);
void InitialSolutionMultires(TPZFMatrix &InitialSol, TPZCompMesh * cmesh);

double ComputeTimeStep(double CFL, int maxLevel, int meshlevel, TPZGeoMesh * gmesh);

TPZGeoMesh * CreateCoarseMesh(int nLevel);
