//$Id: main.cpp,v 1.8 2006-03-16 01:56:04 tiago Exp $

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
#include "TPZTimer.h"

using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;
using namespace std;

#include <set>

int main(){
   
//   TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);
  TPZShapeLinear::fOrthogonal = TPZShapeLinear::Chebyshev;
//  TPZShapeLinear::fOrthogonal = TPZShapeLinear::Jacobi;  
//  TPZShapeLinear::fOrthogonal = TPZShapeLinear::Legendre;  

  cout << "h p" << endl;
  int h, p;
  cin >> h;
  cin >> p;
  cout << "h = " << h << " - p = " << p << endl;  

  TPZCompMesh * cmesh = BarraTracionadaGirada/*BarraTracionadaNeumann*/(h, p);
  TPZGeoMesh *gmesh = cmesh->Reference();
 
  TPZAnalysis an(cmesh);

//#define direct  
#ifdef direct
  /*TPZFrontStructMatrix <TPZFrontSym> */TPZSkylineStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ECholesky);
  an.SetSolver(step);
#endif
  
#define iter;
#ifdef iter
  cout << "ITER_SOLVER" << endl;  
  /*TPZFStructMatrix*/ TPZSkylineStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZStepSolver step( full.Create() );
  an.SetSolver(step);  
  REAL tol = 1.e-8;

#define Precond
#ifdef Precond 

  #define SSORPRECOND
  #ifdef SSORPRECOND
     TPZStepSolver precond( full.Create() );
     REAL precondtol = 1.E-20;
//     precond.SetSSOR( 1, 1., precondtol, 0);     
     precond.SetJacobi(1, precondtol, 0);
     step.ShareMatrix( precond );  
     step.SetCG( 20000, precond, tol, 0 );
//     step.SetGMRES( 2000, 20, precond, tol, 0);     
     cout << "SSOR/JACOBI PRECOND" << endl;
  #else

      TPZMatrixSolver * precond = an.BuildPreconditioner(TPZAnalysis::EBlockJacobi , false);
      step.SetCG( 2000, *precond, tol, 0 );
//      step.SetGMRES( 2000, 20, *precond, tol, 0);
      delete precond;
  #endif

#else
// Sem pre-condicionador 
     TPZCopySolve precond( full.Create() );
     step.ShareMatrix( precond );  
     step.SetCG( 2000, precond, tol, 0 );
//     step.SetGMRES( 2000, 20, precond, tol, 0);     
//     step.SetSSOR( 18000, 1., tol, 0);
     cout << "SEM PRECOND" << endl;
#endif

     an.SetSolver(step);
#endif  
  
  an.Run();
  std::cout << "Numero de equacoes = " << an.Mesh()->NEquations() << std::endl;
  cout << "h = " << h << " - p = " << p << endl;  
  std::cout.flush();
  an.Run();

/**** Aqui faz DX ****/
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(1);
  scalnames[0] = "Stress1";
  vecnames[0] = "PrincipalStress";
  std::stringstream filedx;
  filedx << "result.dx";
  an.DefineGraphMesh(3,scalnames,vecnames,&(filedx.str()[0]));
  an.PostProcess(1);
  
  delete cmesh;
  delete gmesh;
  return 0;
}

