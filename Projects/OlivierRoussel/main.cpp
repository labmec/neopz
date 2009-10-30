//$Id: main.cpp,v 1.7 2009-10-30 13:18:40 cesar Exp $

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
#include "adapt.h"

using namespace std;

int main(int argc, char *argv[])
{
  if (argc > 1) 
  {
	std::string logpath ( argv[1] );
    cout << "initializing LOG usign the following configuration file " << logpath << endl;
    InitializePZLOG ( logpath );	
  }
  else
  {
    cout << "initializing LOG\n";
	InitializePZLOG();
  }
  const int L = 5;

  REAL timeStep;
//   TPZCompMesh * cmesh = CreateMeshLaxAndSod(L,timeStep);
//   TPZCompMesh * cmesh = CreateMeshLax2D(L,timeStep);

  cout << "generating interpolation space...\n";
  TPZCompMesh * cmesh = CreateMeshMultires(gmesh);

  
//  cout << "loading dummy solution...\n";
  
//  LoadDummySolution ( cmesh );
  
 // cout << "entering adaptive procedure\n";
  //GetAdaptedMesh ( cmesh );
  
 // cmesh->Print();
 // cout << "finish!\n";
 // return EXIT_SUCCESS;
  

  timeStep = ComputeTimeStep(0.01,L,L,cmesh->Reference());


#ifdef DEBUG
{
  ofstream malhas("malhas.txt");
  cmesh->Reference()->Print(malhas);
  cmesh->Print(malhas);
}
#endif

  TPZExplFinVolAnal an(cmesh, cout);

  InitializeSolver(an);
  const double PhysicalTime = 0.05;
  int niter = PhysicalTime/timeStep+1;
  cout << "\nnequations = " << cmesh->NEquations();
  cout << "\nNiter = " << niter << "\n";

  TPZFMatrix InitialSol;
//   InitialSolutionLaxAndSod(InitialSol,cmesh);
//   InitialSolutionLax2D(InitialSol,cmesh);
  InitialSolutionLinearConvection(InitialSol,cmesh);
//   InitialSolutionMultires(InitialSol,cmesh);
  an.SetInitialSolution(InitialSol);

  an.Set(timeStep,niter,1e-10);

//  an.SetSaveFrequency(niter/5,0);
//  TPZVec<string> scal(3-2),vec(0);
//  scal[0] = "density";
//  scal[1] = "energy";
//  scal[2] = "Mach";
//  stringstream nome; nome << "testeL"<<".dx";
//  an.DefineGraphMesh(3,scal,vec,nome.str());
//   an.Run();	

	double Epsl = 0.125 * 0.1*5.0;
  an.MultiResolution( Epsl );


  delete cmesh;
  delete gmesh;
  return EXIT_SUCCESS;  
}

