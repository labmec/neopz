//$Id: main.cpp,v 1.2 2005-12-13 19:09:08 tiago Exp $

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

#include <time.h>
#include <stdio.h>
#include "meshes.h"

using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;
using namespace std;

int main(){
   
//   TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);
  TPZCompMesh * cmesh = /*VigaEngastada*/BarraTracionada(1, 2);
  TPZGeoMesh *gmesh = cmesh->Reference();
 
  TPZAnalysis an(cmesh);

  /*TPZParFrontStructMatrix <TPZFrontSym>*/ TPZFStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ECholesky/*ELU*/);
  an.SetSolver(step);

  an.Run();

/**** Aqui faz DX ****/
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(1);
  scalnames[0] = "DisplacementX";
  vecnames[0] = "Displacement";
  std::stringstream filedx;
  filedx << "result.dx";
  an.DefineGraphMesh(3,scalnames,vecnames,&(filedx.str()[0]));
  an.PostProcess(0);
  
  delete cmesh;
  delete gmesh;
  return 0;
}

