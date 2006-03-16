//$Id: meshes.cpp,v 1.5 2006-03-16 01:53:09 tiago Exp $

#include "meshes.h"

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

#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
#include "pzintel.h"
#include "pzbndcond.h"
#include "pzpoisson3d.h"
#include "pzcoupledtransportdarcy.h"
#include "pzcoupledtransportdarcyBC.h"
#include "TPZShapeDisc.h"
#include "pzgmesh.h"
#include "pzshapelinear.h"

using namespace pzshape;

void SetPOrder(int p){
  TPZCompEl::gOrder = p;
}

TPZCompMesh * CreateSimpleMeshWithExactSolution(int h, int p){

  TPZCompEl::gOrder = p;

  REAL co[4][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}};
  int indices[1][4] = {{0,1,2,3}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 4;
  int nelem = 1;
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
  
  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],5,-1,*gmesh); // right
  TPZGeoElBC gbc3(elvec[0],6,-1,*gmesh); // right
  TPZGeoElBC gbc4(elvec[0],7,-1,*gmesh); // top
  
  //Refinamento uniforme
  for(int ref = 0; ref < h; ref++){
    TPZManVector<TPZGeoEl *> filhos;
    int n = gmesh->NElements();    
    for(int i = 0; i < n; i++){
      TPZGeoEl * gel = gmesh->ElementVec()[i];
      if (gel->Dimension() == 2) gel->Divide(filhos);
    }//for i
  }//ref

  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZCoupledTransportDarcy * mat = new TPZCoupledTransportDarcy(1, 2, 3, 2);
  TPZMatPoisson3d * FirstEq  = mat->GetMaterial(0);
  TPZMatPoisson3d * SecondEq = mat->GetMaterial(1);
  FirstEq->SetInternalFlux(0.);
  SecondEq->SetForcingFunction(Forcing2);

  TPZManVector<REAL,2> convdir(2,0.);
  FirstEq->SetParameters(1., 0., convdir);
  SecondEq->SetParameters(1./2., 0., convdir);
  mat->SetAlpha(1./3.);
  
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc;  
  val2.Zero();
  TPZCoupledTransportDarcyBC * BC = mat->CreateBC(-1);
  
  bc = FirstEq->CreateBC(-2, 0,val1,val2);
  bc->SetForcingFunction(Dirichlet1);
  BC->SetMaterial(0, bc);
  
  val2(0,0) = 1.;
  bc = SecondEq->CreateBC(-3, 0, val1, val2);
  BC->SetMaterial(1, bc);

  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(BC);

  cmesh->SetAllCreateFunctionsContinuous();  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();
  
  return cmesh;
}//end of method

void ExactSol_p(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv){
  REAL x = pt[0];
  REAL y = pt[1];
  sol.Resize(1);
  sol[0] = x*y;
  deriv(0,0) = y;
  deriv(1,0) = x;  
}

void ExactSol_u(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv){
  REAL x = pt[0];
  REAL y = pt[1];
  sol.Resize(1);
  sol[0] = x*y*(x-1.)*(y-1.)+1.;
  deriv.Resize(2,1);
  deriv(0,0) = (-1.+x)*(-1.+y)*y+x*(-1.+y)*y;
  deriv(1,0) = (-1.+x)*x*(-1.+y)+(-1.+x)*x*y;
}

void Dirichlet1(TPZVec<REAL> &pt, TPZVec<REAL> &force){
  REAL x = pt[0];
  REAL y = pt[1];
  force.Resize(1);
  force[0] = x*y;
}
  
void Forcing2(TPZVec<REAL> &pt, TPZVec<REAL> &force){
  REAL x = pt[0];
  REAL y = pt[1];
  force.Resize(1);
  REAL val = 1./2.*(-2.*(-1.+x)*x-2.*(-1.+y)*y)
            +1./3.*(-x*((-1.+x)*x*(-1.+y)+(-1.+x)*x*y)-y*((-1.+x)*(-1.+y)*y+x*(-1.+y)*y));
  force[0] = -1. * val;
}

TPZCompMesh * CheckBetaNonConstant(int h, int p){
  
  TPZCompEl::gOrder = p;

  REAL co[4][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}};
  int indices[1][4] = {{0,1,2,3}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 4;
  int nelem = 1;
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
  
  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],5,-1,*gmesh); // right
  TPZGeoElBC gbc3(elvec[0],6,-1,*gmesh); // right
  TPZGeoElBC gbc4(elvec[0],7,-1,*gmesh); // top
  
  //Refinamento uniforme
  for(int ref = 0; ref < h; ref++){
    TPZManVector<TPZGeoEl *> filhos;
    int n = gmesh->NElements();    
    for(int i = 0; i < n; i++){
      TPZGeoEl * gel = gmesh->ElementVec()[i];
      if (gel->Dimension() == 2) gel->Divide(filhos);
    }//for i
  }//ref

  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d * mat = new TPZMatPoisson3d(1,2);

  TPZManVector<REAL,2> convdir(2,0.);
  mat->SetParameters(1., 0., convdir);
  mat->SetForcingFunction(Force);
  
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc;  
  val2.Zero();
  TPZBndCond * BC = mat->CreateBC(-1, 0, val1, val2);
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(BC);

  cmesh->SetAllCreateFunctionsContinuous();  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
//   cmesh->CleanUpUnconnectedNodes();
//   cmesh->ExpandSolution();
  
  return cmesh;
}
void SolExata(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv){
    REAL x = pt[0]; REAL y = pt[1];
    const REAL Pi = 4. * atan(1.);
    sol.Resize(1);
    sol[0] = sin(Pi*x)*sin(Pi*y);
    deriv.Resize(2,1);
    deriv(0,0) = Pi * cos(Pi*x) * sin(Pi*y);
    deriv(1,0) = Pi * cos(Pi*y) * sin(Pi*x);
}

void Force(TPZVec<REAL> &pt, TPZVec<REAL> &force){
  REAL x = pt[0];
  REAL y = pt[1];
  force.Resize(1);
  const REAL Pi = 4. * atan(1.);
  REAL val = Pi*x*cos(Pi*y)*sin(Pi*x)-Pi*y*cos(Pi*x)*sin(Pi*y)+2.*Pi*Pi*sin(Pi*x)*sin(Pi*y);
  force[0] = -1. * val;
}


TPZCompMesh * CreateSimpleMeshWithExactSolution2(int h, int p){

  TPZCompEl::gOrder = p;

  REAL co[4][2] = {{0.,0.},{1.,0.},{1.,1.},{0.,1.}};
  int indices[1][4] = {{0,1,2,3}};
  TPZGeoEl *elvec[1];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 4;
  int nelem = 1;
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
  
  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],5,-1,*gmesh); // right
  TPZGeoElBC gbc3(elvec[0],6,-1,*gmesh); // right
  TPZGeoElBC gbc4(elvec[0],7,-1,*gmesh); // top
  
  //Refinamento uniforme
  for(int ref = 0; ref < h; ref++){
    TPZManVector<TPZGeoEl *> filhos;
    int n = gmesh->NElements();    
    for(int i = 0; i < n; i++){
      TPZGeoEl * gel = gmesh->ElementVec()[i];
      if (gel->Dimension() == 2) gel->Divide(filhos);
    }//for i
  }//ref

  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZCoupledTransportDarcy * mat = new TPZCoupledTransportDarcy(1, 2, 3, 2);
  TPZMatPoisson3d * FirstEq  = mat->GetMaterial(0);
  TPZMatPoisson3d * SecondEq = mat->GetMaterial(1);
  FirstEq->SetInternalFlux(0.);
  SecondEq->SetForcingFunction(Forcing22_eMenos6/*Forcing22*/);

  TPZManVector<REAL,2> convdir(2,0.);
  FirstEq->SetParameters(1., 0., convdir);
  SecondEq->SetParameters(1.e-06/*1./2.*/, 0., convdir);
  mat->SetAlpha(1./3.);
  
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc;  
  val2.Zero();
  TPZCoupledTransportDarcyBC * BC = mat->CreateBC(-1);
  
  bc = FirstEq->CreateBC(-2, 0,val1,val2);
  bc->SetForcingFunction(Dirichlet1);
  BC->SetMaterial(0, bc);
  
  val2(0,0) = 1.;
  bc = SecondEq->CreateBC(-3, 0, val1, val2);
  BC->SetMaterial(1, bc);

  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(BC);

  cmesh->SetAllCreateFunctionsDiscontinuous();  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();
  
  return cmesh;
}//end of method

void ExactSol_u2(TPZVec<REAL> &pt, TPZVec<REAL> &sol, TPZFMatrix &deriv){
  REAL x = pt[0];
  REAL y = pt[1];
  sol.Resize(1);
  sol[0] = exp(x)*x*y*(x-1.)*(y-1.)+1.;
  deriv.Resize(2,1);
  deriv(0,0) = exp(x)*(-1.+x+x*x)*(-1.+y)*y;
  deriv(1,0) = exp(x)*(-1.+x)*x*(-1.+2.*y);
}

void Forcing22(TPZVec<REAL> &pt, TPZVec<REAL> &force){
  REAL x = pt[0];
  REAL y = pt[1];
  force.Resize(1);
  REAL val = -1./6.*exp(x)*(
             -2.*(-1.+y)*y*y+x*x*x*(-2.+4.*y)
             +x*x*(8.-7.*y+y*y+2.*y*y*y)
             +x*(-6.-9.*y+7.*y*y+2.*y*y*y)
             );
  force[0] = -1. * val;
}

void Forcing22_eMenos6(TPZVec<REAL> &pt, TPZVec<REAL> &force){
  REAL x = pt[0];
  REAL y = pt[1];
  force.Resize(1);
  REAL val = -exp(x)*x*(-2.-3.*y+3.*y*y+x*(2.-y+y*y))*1.e-06
             -1./3.*exp(x)*(-(-1.+y)*y*y+x*(-1.+y)*y*y+x*x*x*(-1.+2.*y)+x*x*(1.-2.*y-y*y+y*y*y));
  force[0] = -1. * val;
}
