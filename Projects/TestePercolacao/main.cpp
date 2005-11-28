//$Id: main.cpp,v 1.2 2005-11-28 17:16:42 tiago Exp $

/**
 * Percolation of water from the fracture into the porous media.
 * September 16, 2005
 */

#include <sstream>
#include "pzvec.h"
#include "pzcmesh.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"
//#include "pzerror.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

//#include "pzintel.h"
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
#include "pzonedref.h"


#include "pzbndcond.h"
#include "pzpoisson3d.h"

#include "pzvisualmatrix.h"

#include "TPZGeoElement.h"
#include "pzgeoel.h"

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

#include <time.h>
#include <stdio.h>
#include "meshes.h"
#include "postprocess.h"

using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;
using namespace std;
 
int main(){   

  TPZCompMesh * cmesh = CreateMesh_TriangularDomain_ComoPhilippeQuer(3, 30, 6);// CreateMeshSingularity(1, 2,  2);
  std::cout << "Numero de elementos = " << cmesh->ElementVec().NElements() << std::endl;
  std::cout << "Numero de equacoes  = " << cmesh->NEquations() << std::endl;
  TPZGeoMesh *gmesh = cmesh->Reference();
 
  TPZAnalysis an(cmesh);

//#define ITER_SOLVER
#ifdef ITER_SOLVER
  cout << "ITER_SOLVER" << endl;  
  TPZSpStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZCopySolve precond( full.Create()  );
  TPZStepSolver step;
  step.ShareMatrix( precond );  
//  step.SetGMRES( 2000000, 30, precond, 1.e-13, 0 ); 
  step.SetCG( 2000000, precond, 1e-15, 0);
  an.SetSolver(step);
#endif

#define FRONTAL
#ifdef FRONTAL  
  TPZParFrontStructMatrix <TPZFrontSym> /* TPZFStructMatrix */ full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ECholesky/*ELU*/);
  an.SetSolver(step);
#endif
  an.Run();
  
   std::ofstream FluxFile("resultado/flux.txt");
   ComputeNormalFlux(FluxFile, *cmesh, -1);
   std::ofstream SolFile("resultado/solution.txt");
   ComputeSolution(SolFile, *cmesh, -1);

std::cout << "\nFEITO\n";
std::cout.flush();
/**** Aqui faz DX ****/
  const int Resolution = 4;
  {
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(0);
  scalnames[0] = "NormKDu";//"KDuDy";//"Solution";
//  vecnames[0] = "Derivative";
  std::stringstream filedx; 
  filedx << "resultado/resultFEM_";
  filedx << "NormKDu.dx";
  an.DefineGraphMesh(2,scalnames,vecnames,&(filedx.str()[0]));
  an.PostProcess(Resolution);
  }
  {
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(0);
  scalnames[0] = "Solution";
//  vecnames[0] = "Derivative";
  std::stringstream filedx;
  filedx << "resultado/resultFEM_";
  filedx << "Solution.dx";
  an.DefineGraphMesh(2,scalnames,vecnames,&(filedx.str()[0]));
  an.PostProcess(Resolution);
  }
  {
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(0);
  scalnames[0] = "KDuDy";
//  vecnames[0] = "Derivative";
  std::stringstream filedx;
  filedx << "resultado/resultFEM_";
  filedx << "KDuDy.dx";
  an.DefineGraphMesh(2,scalnames,vecnames,&(filedx.str()[0]));
  an.PostProcess(Resolution);  
  }
  {
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(0);
  scalnames[0] = "KDuDx";
//  vecnames[0] = "Derivative";
  std::stringstream filedx;
  filedx << "resultado/resultFEM_";
  filedx << "KDuDx.dx";
  an.DefineGraphMesh(2,scalnames,vecnames,&(filedx.str()[0]));
  an.PostProcess(Resolution);  
  }  
  delete cmesh;
  delete gmesh;
  return 0; 
}

