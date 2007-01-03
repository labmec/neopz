#include "pzvec.h"

#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzgnode.h"

#include "pzmatrix.h"

#include "pzelasmat.h"
#include "pzbndcond.h"
//#include "pzplaca.h"
#include "pzmat2dlin.h"


const int nstate = 1;
const REAL onethird = 0.33333333333333333;
const REAL PI = 3.141592654;

static void Neumann2(TPZVec<REAL> &x, TPZVec<REAL> &force);
static void Neumann3(TPZVec<REAL> &x, TPZVec<REAL> &force);
static void Neumann4(TPZVec<REAL> &x, TPZVec<REAL> &force);
static void Neumann5(TPZVec<REAL> &x, TPZVec<REAL> &force);

void Neumann2(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  	REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
  	REAL theta = atan2(x[1],x[0]);
  	REAL rexp = pow(r,onethird);
  	f[0] = -onethird*cos(onethird*(PI/2.-2.*theta))/(rexp*rexp);
}
void Neumann3(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  	REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
  	REAL theta = atan2(x[1],x[0]);
  	REAL rexp = pow(r,onethird);
  	f[0] = onethird*sin(onethird*(PI/2.-2.*theta))/(rexp*rexp); 
}
void Neumann4(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  	REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
  	REAL theta = atan2(x[1],x[0]);
  	REAL rexp = pow(r,onethird);
  	f[0] = onethird*cos(onethird*(PI/2.-2.*theta))/(rexp*rexp);
}
void Neumann5(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  	REAL r = sqrt(x[0]*x[0]+x[1]*x[1]);
  	REAL theta = atan2(x[1],x[0]);
  	REAL rexp = pow(r,onethird);
  	f[0] = -onethird*sin(onethird*(PI/2.-2.*theta))/(rexp*rexp);
}

static TPZCompMesh *CreateMesh();


//*************************************
//************Option 1*****************
//*******L Shape Quadrilateral*********
//*************************************
TPZCompMesh *CreateMesh() {
  REAL co[8][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.}};
  int indices[3][4] = {{0,1,2,3},{0,3,4,5},{0,5,6,7}};
  TPZGeoEl *elvec[3];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 8;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZVec<REAL> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    TPZGeoNode pznode(nod,coord,*gmesh);
    gmesh->NodeVec()[nodind] = pznode;
  }

  int el;
  int nelem = 3;
  for(el=0; el<nelem; el++) {
    TPZVec<int> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
    //    elvec[el] = new TPZGeoElQ2d(el,nodind,1,*gmesh);
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  }
  
  gmesh->BuildConnectivity2();
  
  TPZVec<TPZGeoEl *> sub;
 // elvec[0]->Divide(sub);
//  elvec[1]->Divide(sub);
//  elvec[2]->Divide(sub);
  
//  TPZGeoElBC gbc;
  
  // bc -1 -> Dirichlet
  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh);
  // bc -2 -> Neumann at the bottom y==-1
  TPZGeoElBC gbc2(elvec[0],5,-2,*gmesh);
  // bc -3 -> Neumann at the right x==1
  TPZGeoElBC gbc3(elvec[0],6,-3,*gmesh);
  
  // bc -3 -> Neumann at the right x==1
  TPZGeoElBC gbc4(elvec[1],5,-3,*gmesh);
  
  // bc -4 -> Neumann at the top y==1
  TPZGeoElBC gbc5(elvec[1],6,-4,*gmesh);
  
  // bc -4 -> Neumann at the top y==1
  TPZGeoElBC gbc6(elvec[2],5,-4,*gmesh);
  
  // bc -5 -> Neumann at the left x==-1
  TPZGeoElBC gbc7(elvec[2],6,-5,*gmesh);
  
  // bc -6 -> Homogeneous Neumann
  TPZGeoElBC gbc8(elvec[2],7,-6,*gmesh);
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  
  TPZAutoPointer<TPZMaterial> mat;
  if(nstate == 2) {
    mat = new TPZElasticityMaterial(1,2.,0.3,1.,1.);
  } else {
    TPZMat2dLin *mat2d = new TPZMat2dLin(1);
    int ist,jst;
    TPZFMatrix xk(nstate,nstate,1.),xc(nstate,nstate,0.),xf(nstate,1,0.);
    for(ist=0; ist<nstate; ist++) {
      if(nstate != 1) xf(ist,0) = 1.;
      for(jst=0; jst<nstate; jst++) {
	if(ist != jst) xk(ist,jst) = 0.;
      }
    }
    mat2d->SetMaterial(xk,xc,xf);
    mat = mat2d;
  }
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZAutoPointer<TPZMaterial>bc[5];
  bc[0] = mat->CreateBC(mat,-1,0,val1,val2);
  int i;
  if(nstate == 1) {
    for(i=1; i<6; i++) {
      bc[i] = mat->CreateBC(mat,-i-1,1,val1,val2);
    }
    bc[1]->SetForcingFunction(Neumann2);
    bc[2]->SetForcingFunction(Neumann3);
    bc[3]->SetForcingFunction(Neumann4);
    bc[4]->SetForcingFunction(Neumann5);
  } else {
    for(i=1; i<6; i++) {
      bc[i] = mat->CreateBC(mat,-i-1,0,val1,val2);
    }
  }
  
  cmesh->InsertMaterialObject(mat);
  for(i=0; i<6; i++) cmesh->InsertMaterialObject(bc[i]);
  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  //  cmesh->ExpandSolution();
  
  return cmesh;
}
