//$Id: main.cpp,v 1.2 2006-01-14 20:03:42 tiago Exp $

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
#include "pzcoupledtransportdarcy.h"

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
  TPZMaterial::gBigNumber= 1.e12;
  TPZCompMesh * cmesh = CheckBetaNonConstant/*CreateSimpleMeshWithExactSolution*/(4,2);
  std::cout << "Numero de elementos = " << cmesh->ElementVec().NElements() << std::endl;
  std::cout << "Numero de equacoes  = " << cmesh->NEquations() << std::endl;
  TPZGeoMesh *gmesh = cmesh->Reference();
 
  TPZAnalysis an(cmesh);

//#define ITER_SOLVER
#ifdef ITER_SOLVER
  cout << "ITER_SOLVER" << endl;  
  TPZSpStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZStepSolver step( full.Create() );
  an.SetSolver(step);  

#define Precond
#ifdef Precond
      TPZMatrixSolver * precond = an.BuildPreconditioner(TPZAnalysis::EElement , false);
      step.SetGMRES( 2000000, 30, *precond, 1.e-15, 0 ); 
      delete precond;
#else
// Sem pre-condicionador 
     TPZCopySolve precond( full.Create() );
     step.ShareMatrix( precond );  
     step.SetGMRES( 2000000, 30, precond, 1.e-15, 0 ); 
     cout << "SEM PRECOND" << endl;
#endif  
  
  an.SetSolver(step);
#endif

#define DIRECT
#ifdef DIRECT  
  /*TPZParFrontStructMatrix*//*TPZFrontStructMatrix <TPZFrontNonSym>*/TPZFStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
   step.SetDirect(ELU);
  an.SetSolver(step);
#endif

  TPZCoupledTransportDarcy::SetCurrentMaterial(0);
  std::cout << "\nCalling an.Run() for FirstEq\n";
  an.Run();
//   an.SetExact(ExactSol_p);
  an.SetExact(SolExata);
  TPZVec<REAL> pos;
  an.PostProcess(pos,std::cout);
  std::cout << "Problem solved\n";
  
  {
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(0);
  scalnames[0] = "Solution";
//  vecnames[0] = "Derivative";
  std::stringstream filedx;
  filedx << "1stEq_";
  filedx << "Solution.dx";
  an.DefineGraphMesh(2,scalnames,vecnames,&(filedx.str()[0]));
  an.PostProcess(4);  
  }      

  return 1;
  
  TPZCoupledTransportDarcy::SetCurrentMaterial(1);
  std::cout << "\nCalling an.Run() for SecondEq\n";
    an.Solution().Zero();
  an.Assemble();
  an.Solution().Zero();
  an.Solve();
  an.SetExact(ExactSol_u);
  an.PostProcess(pos,std::cout);  
  std::cout << "Problem solved\n";  
   
  {
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(0);
  scalnames[0] = "Solution";
//  vecnames[0] = "Derivative";
  std::stringstream filedx;
  filedx << "2ndEq_";
  filedx << "Solution.dx";
  an.DefineGraphMesh(2,scalnames,vecnames,&(filedx.str()[0]));
  an.PostProcess(4);  
  }    
    
  delete cmesh;
  delete gmesh;
  return 0; 
}

