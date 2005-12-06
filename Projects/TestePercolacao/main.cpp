//$Id: main.cpp,v 1.3 2005-12-06 13:37:24 tiago Exp $

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
#include "pzskylstrmatrix.h"

#include <time.h>
#include <stdio.h>
#include "meshes.h"
#include "postprocess.h"

using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;
using namespace std;
 
int main22(){
  
  {
  REAL L = 1.;
  const int n = 100;
  std::ofstream file("KRebocoVal1.txt");
  TPZVec<REAL> x(2, 0.0);
  TPZFMatrix result(1,1,0.);
  file << "{ ";
  for(int i = 0; i <= n; i++){
    if (i == 80){
      std::cout << "Para ai, tiago";
    }
    x[0] = i * (L / n);
    KRebocoVal1(x, result);
    file << "{ " << x[0] << ", " << result(0,0) << " }";    
    if (i != n) file << ", ";
  }
  file << " }; ";
  }
  
  {
  REAL L = 1.;
  const int n = 100;
  std::ofstream file("KRebocoVal2.txt");
  TPZVec<REAL> x(2, 0.0);
  TPZVec<REAL> result(1,0.);
  file << "{ ";
  for(int i = 0; i <= n; i++){
    if (i == 80){
      std::cout << "Para ai, tiago";
    }
    x[0] = i * (L / n);
    KRebocoVal2(x, result);
    file << "{ " << x[0] << ", " << result[0] << " }";    
    if (i != n) file << ", ";
  }
  file << " }; ";
  }  
  
}

int main(){   

  TPZCompMesh * cmesh = CreateMesh_ComoPhilippeQuer_Adimensional_Sem_Simetria(1,1,2);//(5, 30, 4);//(5,25,2);
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
  TPZParFrontStructMatrix <TPZFrontNonSym> /*TPZFStructMatrix*/  /*TPZSkylineStructMatrix*/ full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
   step.SetDirect(/*ECholesky*/ELU);
  an.SetSolver(step);
#endif
  std::cout << "\nCalling an.Run()\n";
  an.Run();
  std::cout << "Problem solved\nStarting post-processing\n";
  
   std::ofstream FluxFile("resultado/flux.txt");
   ComputeNormalFlux(FluxFile, *cmesh, -1);
   std::ofstream SolFile("resultado/solution.txt");
   ComputeSolution(SolFile, *cmesh, -1);

std::cout << "\nFEITO\n";
std::cout.flush();
/**** Aqui faz DX ****/
  const int Resolution = 0;
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

