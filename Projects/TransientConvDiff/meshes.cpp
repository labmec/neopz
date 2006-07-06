//$Id: meshes.cpp,v 1.1 2006-07-06 15:48:25 tiago Exp $

#include "meshes.h"


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

#include "pzadmchunk.h"


#include "pzbndcond.h"
#include "pzpoisson3d.h"

#include "pzvisualmatrix.h"

#include <time.h>
#include <stdio.h>

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
#include "pztransientanalysis.h"

 using namespace pzgeom;
 using namespace pzshape;
 using namespace pzrefine;
 using namespace std;

const double eps = 1.0;

///////////////////////////////////FUNCTIONS///////////////////////////////////////////

void ForcingFunction(TPZVec<REAL> &pto, TPZVec<REAL> &force){
  REAL T = TPZTransientAnalysis::gTime;
  REAL x = pto[0];
  REAL y = pto[1];  
  force.Resize(1);  
  REAL val = exp(-T)*(-(-1.+T)*(-1.+x)*x*(-1.+y)*y-2.*eps*T*(-x+x*x+(-1.+y)*y));  
  force[0] = val;
  //*-1 pq o pzpoisson3d implementa -Epsilon Laplac(u) + div(beta*u) = -force  
  force[0] *= -1.0;
}

void ExactSolution(TPZVec<REAL> &pto, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  REAL T = TPZTransientAnalysis::gTime;
  REAL x = pto[0];
  REAL y = pto[1];
  u.Resize(1);
  deriv.Resize(2,1);
  
  REAL solucao = x*(x-1.)*y*(y-1.)*T*exp(-T);
  u[0] = solucao;
  

  REAL dx = exp(-T)*T*(-1.+x)*(-1.+y)*y+exp(-T)*T*(-1.+y)*x*y;
  REAL dy = exp(-T)*T*(-1.+x)*(-1.+y)*x+exp(-T)*T*(-1.+x)*x*y;
  
  deriv(0,0) = dx;
  deriv(1,0) = dy;

//  cout << pto << "\t" << u << "\t" << deriv << endl;

}

void ForcingFunction2(TPZVec<REAL> &pto, TPZVec<REAL> &force){
  REAL T = TPZTransientAnalysis::gTime;
  REAL x = pto[0];
  REAL y = pto[1];  
  force.Resize(1);  
  REAL val = (-1.+x)*x*(-1.+y)*y-2.*eps*T*((-1.+x)*x+(-1.+y)*y);
  force[0] = val;
  //*-1 pq o pzpoisson3d implementa -Epsilon Laplac(u) + div(beta*u) = -force  
  force[0] *= -1.0;
}

void ExactSolution2(TPZVec<REAL> &pto, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  REAL T = TPZTransientAnalysis::gTime;
  REAL x = pto[0];
  REAL y = pto[1];
  u.Resize(1);
  deriv.Resize(2,1);
  
  REAL solucao = x*(x-1.)*y*(y-1.)*T;
  u[0] = solucao;
  

  REAL dx = T*(-1.+x)*(-1.+y)*y+T*(-1.+y)*x*y;
  REAL dy = T*(-1.+x)*(-1.+y)*x+T*(-1.+x)*x*y;
  
  deriv(0,0) = dx;
  deriv(1,0) = dy;

//  cout << pto << "\t" << u << "\t" << deriv << endl;

}

#include "pztransientmat.h"
#include "pznonlinearpoisson3d.h"
TPZCompMesh *CreateMesh(int h) {

  REAL co[9][2] = {{0.5,0.5},{0.5,0.},{1.,0.},{1.,0.5},{1.,1.},{0.5,1.},{0.,1.},{0.,0.5},{0.,0.}};
  int indices[4][4] = {{1,2,3,0},{0,3,4,5},{7,0,5,6},{8,1,0,7}};
  TPZGeoEl *elvec[4];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 9;
  int nelem = 4;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<nelem; el++) {
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  }

  gmesh->BuildConnectivity();

  TPZVec<TPZGeoEl*> filhos;
  for(int i = 0; i < h; i++){
     int n = gmesh->NElements();
     for(int j = 0; j < n; j++){
	if (gmesh->ElementVec()[j]->Dimension() == 2) gmesh->ElementVec()[j]->Divide(filhos);
     }
  }

  //TPZGeoElBC gbc;

  TPZGeoElBC gbc1(elvec[0],4,-3,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],5,-3,*gmesh); // right
  TPZGeoElBC gbc3(elvec[1],5,-3,*gmesh); // right
  TPZGeoElBC gbc4(elvec[1],6,-3,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],6,-3,*gmesh); // top
  TPZGeoElBC gbc6(elvec[2],7,-3,*gmesh); // left
  TPZGeoElBC gbc7(elvec[3],7,-3,*gmesh); // left
  TPZGeoElBC gbc8(elvec[3],4,-3,*gmesh); // bottom
    
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new /*TPZMatPoisson3d(1,2);*/TPZTransientMaterial< TPZNonLinearPoisson3d >(1, 2, 0.2); dynamic_cast<TPZNonLinearPoisson3d*>(mat)->SetReferred(false);
  mat->SetForcingFunction(ForcingFunction);
  mat->SetNonSymmetric();
  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  REAL beta = 0.0;
  mat->SetParameters(eps, beta, convdir);
  int nstate = 1;
  
  
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc;

  //Dirichlet nulo
  val2.Zero();
  bc = mat->CreateBC(-3, 0,val1,val2);
 
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bc);

  cmesh->SetAllCreateFunctionsContinuous();
//  cmesh->SetAllCreateFunctionsDiscontinuous();
  
  cmesh->AutoBuild();

  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  //  cmesh->ExpandSolution();
  
  return cmesh;
}

