//$Id: main.cpp,v 1.4 2006-01-09 13:56:10 tiago Exp $

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
  
  TPZCompMesh * cmesh = BarraTracionadaGirada(2, 2);
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
  TPZSkylineStructMatrix /*TPZFStructMatrix*/ full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZStepSolver step( full.Create() );
  an.SetSolver(step);  
  REAL tol = 1.e-8;

#define Precond
#ifdef Precond
      TPZMatrixSolver * precond = an.BuildPreconditioner(TPZAnalysis::EElement , false);
      step.SetCG( 2000, *precond, tol, 0 );
      delete precond;
#else
// Sem pre-condicionador 
     TPZCopySolve precond( full.Create() );
     step.ShareMatrix( precond );  
     step.SetCG( 2000, precond, tol, 0 );
     cout << "SEM PRECOND" << endl;
#endif

     an.SetSolver(step);
#endif  
  
  an.Run();
  std::cout << "Numero de equacoes = " << an.Mesh()->NEquations() << std::endl;
  std::cout.flush();
   
  delete cmesh;
  delete gmesh;
  return 0;
}

