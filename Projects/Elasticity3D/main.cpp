//$Id: main.cpp,v 1.6 2006-01-14 20:02:07 tiago Exp $

/**
 * Validation test of TPZElasticity3D material
 * 24/11/2003
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

using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;
using namespace std;

#include <set>
int main22(){

  TPZFMatrix mat(3,3);
  mat(0,0) = 1999.; mat(0,1) = 2.; mat(0,2) = 3.;
  mat(1,0) = 2.; mat(1,1) = 4.; mat(1,2) = 5.;
  mat(2,0) = 3.; mat(2,1) = 5.; mat(2,2) = 6.;
  mat.Print("Matriz original", cout, EMathematicaInput);

/*  int nn;
  cin >> nn;  
  const int n = nn;//50;
  TPZFMatrix mat(n,n);
  for(int i = 0; i < n; i++)
    for(int j = 0; j <= i; j++){
      mat(i,j) = exp(i*j*0.0001)*(i+1)*(j+1);
      mat(j,i) = exp(i*j*0.0001)*(i+1)*(j+1);
    }    
   mat.Print("Matriz original", cout, EMathematicaInput);   */   
  

  
  int niter = 100000;
  REAL tol = 1.e-200;
  TPZVec< REAL > Sort;
//   mat.SolveEigenvaluesJacobi(niter, tol, &Sort);
//   cout << "Sort:\n";
//   cout << "NElements = " << Sort.NElements() << endl;
//   for(int i = 0; i < Sort.NElements(); i++) cout << Sort[i] << endl;
// 
//   cout << "\nniter = " << niter << " - res = " << tol << endl;
//   mat.Print("Matriz decomposta", cout, EMathematicaInput);
  
  
  TPZFMatrix EigenVec;
  mat.SolveEigensystemJacobi(niter, tol, Sort, EigenVec);
  cout << "AutoValores:" << endl;
  for(int i = 0; i < Sort.NElements(); i++) cout << Sort[i] << endl;
  
  cout << "AutoVetores:" << endl;
  for(int i = 0; i < EigenVec.Rows(); i++){
    for(int j = 0; j < EigenVec.Cols(); j++) cout << EigenVec(i,j) << " , ";
    cout << " fim" << endl;
  }
  
  
}

int main(){
   
  TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);
  TPZShapeLinear::fOrthogonal = TPZShapeLinear::Jacobi;
  
  TPZCompMesh * cmesh = /*BarraTracionadaNeumann*/BarraTracionada/*Girada*//*VigaEngastada*//*VigaEngastadaForcaVolume*/(1, 3);
  
  TPZGeoMesh *gmesh = cmesh->Reference();
 
  TPZAnalysis an(cmesh);

#define direct  
#ifdef direct
  TPZFrontStructMatrix <TPZFrontSym> /*TPZSkylineStructMatrix*/ /*TPZFStructMatrix*/ full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ECholesky/*ELU*/);
  an.SetSolver(step);
#endif
  
//#define iter;
#ifdef iter
  cout << "ITER_SOLVER" << endl;  
  /*TPZFStructMatrix*/ TPZSkylineStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZStepSolver step( full.Create() );
  an.SetSolver(step);  
  REAL tol = 1.e-14;
  TPZStepSolver precond( full.Create() );
  REAL precondtol = 1.E-20;
  precond.SetSSOR( 1, 1., precondtol, 0);     
//  precond.SetJacobi(1, precondtol, 0);
  step.ShareMatrix( precond );  
  step.SetCG( 20000, precond, tol, 0 );
  cout << "SSOR/JACOBI PRECOND" << endl;


     an.SetSolver(step);
#endif  
  
  std::cout << "Numero de equacoes = " << an.Mesh()->NEquations() << std::endl;
  std::cout.flush();
  an.Run();

/**** Aqui faz DX ****/
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(1);
  scalnames[0] = "Strain1";
  vecnames[0] = "PrincipalStrain";
  std::stringstream filedx;
  filedx << "result.dx";
  an.DefineGraphMesh(3,scalnames,vecnames,&(filedx.str()[0]));
  an.PostProcess(0);
  
  delete cmesh;
  delete gmesh;
  return 0;
}

