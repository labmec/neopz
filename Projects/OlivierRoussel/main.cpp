//$Id: main.cpp,v 1.5 2009-10-09 15:14:19 fortiago Exp $

#include "malhas.h"
#include "MultiResMesh.h"
#include "pzlog.h"
#include "pzmatrix.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "pzgeoel.h"
#include <sstream>
#include "TPZGeoElement.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzvec.h"
#include "pzcmesh.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "pzmatrix.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "TPZCopySolve.h"
#include "TPZStackEqnStorage.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"    		
#include "pzbndcond.h"
#include "pzpoisson3d.h"
#include "pzvisualmatrix.h"
#include "TPZGeoElement.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "tpzchangeel.h"
#include "TPZInterfaceEl.h"
#include <time.h>
#include <stdio.h>
#include <set>
#include "tools.h"
#include <TPZTimer.h>
#include "pzexplfinvolanal.h"

using namespace std;

int main(){

  InitializePZLOG();
  const int L = 4;
  REAL timeStep;
//   TPZCompMesh * cmesh = CreateMeshLaxAndSod(L,timeStep);
//   TPZCompMesh * cmesh = CreateMeshLax2D(L,timeStep);
//   TPZCompMesh * cmesh = CreateMeshLinearConvection(L,timeStep);
  TPZGeoMesh * gmesh = CreateCoarseMesh(L);
  TPZCompMesh * cmesh = CreateMeshMultires(gmesh);
  timeStep = ComputeTimeStep(0.5,L,L,gmesh);

#ifdef DEBUG
{
  ofstream malhas("malhas.txt");
  cmesh->Reference()->Print(malhas);
  cmesh->Print(malhas);
}
#endif

  TPZExplFinVolAnal an(cmesh, cout);

  InitializeSolver(an);
  const double PhysicalTime = 0.5;
  int niter = PhysicalTime/timeStep+1;
  cout << "\nnequations = " << cmesh->NEquations();
  cout << "\nNiter = " << niter << "\n";

  TPZFMatrix InitialSol;
//   InitialSolutionLaxAndSod(InitialSol,cmesh);
//   InitialSolutionLax2D(InitialSol,cmesh);
//   InitialSolutionLinearConvection(InitialSol,cmesh);
  InitialSolutionMultires(InitialSol,cmesh);
  an.SetInitialSolution(InitialSol);

  an.Set(timeStep,niter,1e-10);
  an.SetSaveFrequency(niter/1,0);
  TPZVec<string> scal(3-2),vec(0);
  scal[0] = "density";
//  scal[1] = "energy";
//  scal[2] = "Mach";
  stringstream nome; nome << "testeL"/*<<L*/<<".dx";
  an.DefineGraphMesh(3,scal,vec,nome.str());
  an.MultiResolution();
//   an.Run();

  delete cmesh;
  delete gmesh;
  return EXIT_SUCCESS;  
}

