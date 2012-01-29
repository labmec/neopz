//$Id: malhas.cpp,v 1.4 2009-11-04 14:13:24 fortiago Exp $

#include "malhas.h"
#include "TPZFakeFunction.h"
#include "tpzautopointer.h"
#include "pzfunction.h"
#include "pzeuler.h"
#include "pzmatrix.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzgeoelrefpattern.h"
#include "tpzquadraticline.h"
#include "pzgeoel.h"
#include <sstream>
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
#include "pzvec.h"
#include "pzcmesh.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"
#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "pzmatrix.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
// #include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "TPZCopySolve.h"
#include "TPZStackEqnStorage.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"       
#include "pzbndcond.h"
#include "pzpoisson3d.h"
#include "pzvisualmatrix.h"
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
#include "tpzchangeel.h"
#include "TPZInterfaceEl.h"
#include <time.h>
#include <stdio.h>
#include <set>
#include "pzeuler.h"

using namespace std;

// #define ROTACAO
#ifdef ROTACAO
void Rotacao(TPZFMatrix &R){
  R.Redim(3,3);
  double Pi = 4.*atan(1.);
  double t = Pi/6.;
  R(0,0) = +cos(t); R(0,1) = +sin(t); R(0,2) = 0.;
  R(1,0) = -sin(t); R(1,1) = +cos(t); R(1,2) = 0.;
  R(2,0) = 0.;      R(2,1) = 0.;      R(2,2) = 1.;
}

void Rotacao2(TPZFMatrix &R){
  R.Redim(3,3);
  double Pi = 4.*atan(1.);
  double t = Pi/6.;
  R(0,0) = +cos(t);  R(0,1) = 0; R(0,2) = sin(t);
  R(1,0) = 0;        R(1,1) = 1; R(1,2) = 0.;
  R(2,0) = -sin(t);  R(2,1) = 0.; R(2,2) = cos(t);
}

// void ApplyRotacao(TPZFMatrix &R, TPZVec<REAL> &coord){
//   double c0 = coord[0];
//   double c1 = coord[1];
//   double c2 = coord[2];
//   coord[0] = c0*R(0,0) + c1*R(0,1) + c2*R(0,2);
//   coord[1] = c0*R(1,0) + c1*R(1,1) + c2*R(1,2);
//   coord[2] = c0*R(2,0) + c1*R(2,1) + c2*R(2,2);
// }
// 
// void ApplyRotacao(TPZVec<REAL> &coord){
//   TPZFNMatrix<9> R;
//   Rotacao(R);
//   ApplyRotacao(R,coord);
// }

void Gira(TPZFMatrix & R, double &u,double &v,double &w){
  double c0 = u;
  double c1 = v;
  double c2 = w;
  u = c0*R(0,0) + c1*R(0,1) + c2*R(0,2);
  v = c0*R(1,0) + c1*R(1,1) + c2*R(1,2);
  w = c0*R(2,0) + c1*R(2,1) + c2*R(2,2);
}

void ApplyRotacao(double &u,double &v,double &w){
  TPZFNMatrix<9> R;
  Rotacao(R);
  Gira(R,u,v,w);
  Rotacao2(R);
  Gira(R,u,v,w);
}

void RotacionarMalha(TPZCompMesh *cmesh){
  ///girando nos
  {
    double v1=0,v2=0,v3=1;
    ApplyRotacao(v1,v2,v3);
    cout << "\nvetor normal do corte do DX: { " << v1 << ", " << v2 << ", " << v3 << " }\n";
    cout.flush();
  }
  TPZGeoMesh *gmesh = cmesh->Reference();
  int nnodes = gmesh->NodeVec().NElements();
  for(int in = 0; in < nnodes; in++){
    double x = gmesh->NodeVec()[in].Coord(0);
    double y = gmesh->NodeVec()[in].Coord(1);
    double z = gmesh->NodeVec()[in].Coord(2);
    ApplyRotacao(x,y,z);
    gmesh->NodeVec()[in].SetCoord(0,x);
    gmesh->NodeVec()[in].SetCoord(1,y);
    gmesh->NodeVec()[in].SetCoord(2,z);
  }
  int nel = cmesh->NElements();
  for(int iel = 0; iel < nel; iel++){
    TPZCompEl *cel = cmesh->ElementVec()[iel];
    if(!cel) continue;
    TPZInterfaceElement *face = dynamic_cast<TPZInterfaceElement*>(cel);
    if(!face) continue;
    face->SetLeftRightElements(face->LeftElementSide(),face->RightElementSide());///para atualizar a normal
  }
}///void
#endif

TPZCompMesh *CreateMeshLaxAndSod(const int L,REAL &timeStep){
  const int n = (1 << L)+1;

  const REAL deltaX = 2./((double)((n-1.)));
  const REAL CFL = 0.3;
  const REAL Base   = 1.;
  const REAL Height = 1.;
  double min = deltaX; if(Base < min) min = Base;  if(Height < min) min = Height;
  timeStep = CFL*min/2.;

  TPZGeoMesh *gmesh = new TPZGeoMesh();
  for(int i = 0; i < n; i++){
    TPZManVector<REAL,3> coord(3);
    int nodind = gmesh->NodeVec().AllocateNewElement();    
    coord[0] = -1.+deltaX*i;    coord[1] = 0.;    coord[2] = 0.;
    gmesh->NodeVec()[nodind].Initialize(nodind,coord,*gmesh);

    nodind = gmesh->NodeVec().AllocateNewElement();    
    coord[0] = -1.+deltaX*i;    coord[1] = Base;    coord[2] = 0.;
    gmesh->NodeVec()[nodind].Initialize(nodind,coord,*gmesh);

    nodind = gmesh->NodeVec().AllocateNewElement();    
    coord[0] = -1.+deltaX*i;    coord[1] = Base;    coord[2] = Height;
    gmesh->NodeVec()[nodind].Initialize(nodind,coord,*gmesh);

    nodind = gmesh->NodeVec().AllocateNewElement();    
    coord[0] = -1.+deltaX*i;    coord[1] = 0.;    coord[2] = Height;
    gmesh->NodeVec()[nodind].Initialize(nodind,coord,*gmesh);
  }

  for(int el = 0; el < n-1; el++){
    TPZManVector<int,8> nodind(8);
    for(int nod=0; nod < 8; nod++) nodind[nod] = nod+4*el;
    int index;
    gmesh->CreateGeoElement(ECube,nodind,1,index);
  }

  TPZManVector<int,4> nodind(4);
  nodind[0] = 0; nodind[1] = 1; nodind[2] = 2; nodind[3] = 3;
  int index;
  gmesh->CreateGeoElement(EQuadrilateral,nodind,-1,index);
  nodind[0] = 4*n-4; nodind[1] = 4*n-3; nodind[2] = 4*n-2; nodind[3] = 4*n-1;
  gmesh->CreateGeoElement(EQuadrilateral,nodind,-1,index);

  gmesh->BuildConnectivity();

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);

  TPZAutoPointer<TPZMaterial> mat = new TPZEulerEquation(1,1.4);
  cmesh->InsertMaterialObject(mat);
  TPZFMatrix val1,val2;
  cmesh->InsertMaterialObject(mat->CreateBC(mat,-1,TPZEulerEquation::EFreeSlip,val1,val2));

  cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetDefaultOrder(0);
  TPZCompElDisc::SetgOrder(0);

  cmesh->AutoBuild();
  TPZCompElDisc::SetTotalOrderShape(cmesh);
  TPZAutoPointer<TPZFunction> fakefunc = new TPZFakeFunction();
  for(int i = 0; i < cmesh->NElements(); i++){
    TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[i]);
    if(disc){
      disc->SetExternalShapeFunction(fakefunc);
    }
  }
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();

  return cmesh;

}

/** initial solution */
void InitialSolutionLaxAndSod(TPZFMatrix &InitialSol, TPZCompMesh * cmesh){
  InitialSol.Redim(cmesh->NEquations(),1);
  InitialSol.Zero();
  for(int iel = 0; iel < cmesh->NElements(); iel++){
    TPZCompEl * cel = cmesh->ElementVec()[iel];
    if(!cel) continue;
    TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cel);
    if(!disc) continue;
    if(disc->NConnects() == 0) continue;
    int bl = disc->Connect(0).SequenceNumber();
    int blpos = cmesh->Block().Position(bl);
    int blocksize = cmesh->Block().Size(bl);
    
    TPZGeoEl * gel = cel->Reference();
    TPZVec<REAL> xi(3), x(3);
    gel->CenterPoint(gel->NSides()-1,xi);
    gel->X(xi,x);
    double rho = 0, rhoU = 0, rhoV = 0, rhoW = 0, rhoE = 0;
    if(x[0] <= 0.){
      rho = 0.445;//*/1.;
      rhoU = 0.311;//*/0.;
      rhoE = 8.928;//*/2.5;
    }
    else{
      rho = 0.5;//*/0.125;
      rhoU = 0.;//*/0.;
      rhoE = 1.4275;//*/0.25;
    }
#ifdef ROTACAO
    ApplyRotacao(rhoU,rhoV,rhoW);
#endif
    InitialSol(blpos+blocksize-20+0,0) = rho;
    InitialSol(blpos+blocksize-20+1,0) = rhoU;
    InitialSol(blpos+blocksize-20+2,0) = rhoV;
    InitialSol(blpos+blocksize-20+3,0) = rhoW;
    InitialSol(blpos+blocksize-20+4,0) = rhoE;
  }

#ifdef ROTACAO
  RotacionarMalha(cmesh);
#endif

}

TPZCompMesh *CreateMeshLax2D(int L, REAL &timeStep){
  const int n = (1 << L);

  const REAL deltaX = 1./((double)((n)));
  const REAL CFL = 0.45;
  const REAL Height = 0.1;
  double min = deltaX; if(Height < min) min = Height;
  timeStep = CFL*min/2.;

  TPZGeoMesh *gmesh = new TPZGeoMesh();
  for(int k = 0; k < 2; k++){
    for(int j = 0; j < n+1; j++){
      for(int i = 0; i < n+1; i++){
        double x = deltaX*i;
        double y = deltaX*j;
        double z = k*Height;
        TPZManVector<REAL,3> coord(3);
        int nodind = gmesh->NodeVec().AllocateNewElement();
        coord[0] = x;    coord[1] = y;    coord[2] = z;
        gmesh->NodeVec()[nodind].Initialize(nodind,coord,*gmesh);
      }///i
    }///j
  }///k

  for(int j = 0; j < n; j++){
    for(int i = 0; i < n; i++){
      int incid[8] = { i+(n+1)*j, (i+1)+(n+1)*j, (i+1)+(n+1)*(j+1), i+(n+1)*(j+1),
                       i+(n+1)*j+(n+1)*(n+1), (i+1)+(n+1)*j+(n+1)*(n+1), (i+1)+(n+1)*(j+1)+(n+1)*(n+1), i+(n+1)*(j+1)+(n+1)*(n+1)
                     };
      TPZManVector<int,8> nodind(8);
      for(int nod=0; nod < 8; nod++) nodind[nod] = incid[nod];
      int index;
      gmesh->CreateGeoElement(ECube,nodind,1,index);

      ///face en bas
//       TPZManVector<int,4> nodindFace(4);
//       for(int nod=0; nod < 4; nod++) nodindFace[nod] = incid[nod];
//       gmesh->CreateGeoElement(EQuadrilateral,nodindFace,-1,index);

      ///face en haut
/*      for(int nod=0; nod < 4; nod++) nodindFace[nod] = incid[nod+4];
      gmesh->CreateGeoElement(EQuadrilateral,nodindFace,-1,index);*/
    }///i
  }///j

  ///south
  for(int j = 0; j < n; j++){
    int incid[4] = { j, j+1,
                     j+1+(n+1)*(n+1),j+(n+1)*(n+1)
                    };
    ///face
    int index;
    TPZManVector<int,4> nodindFace(4);
    for(int nod=0; nod < 4; nod++) nodindFace[nod] = incid[nod];
    gmesh->CreateGeoElement(EQuadrilateral,nodindFace,-1,index);

  }///j

  ///north
  for(int j = 0; j < n; j++){
    int incid[4] = { (n+1)*n+j,(n+1)*n+j+1,
                     (n+1)*n+j+1+(n+1)*(n+1), (n+1)*n+j+(n+1)*(n+1) 
                    };
    ///face
    int index;
    TPZManVector<int,4> nodindFace(4);
    for(int nod=0; nod < 4; nod++) nodindFace[nod] = incid[nod];
    gmesh->CreateGeoElement(EQuadrilateral,nodindFace,-1,index);

  }///j

  ///west
  for(int i = 0; i < n; i++){
    int incid[4] = { i*(n+1), (i+1)*(n+1),
                     (i+1)*(n+1)+(n+1)*(n+1), i*(n+1)+(n+1)*(n+1)
                    };
    ///face
    int index;
    TPZManVector<int,4> nodindFace(4);
    for(int nod=0; nod < 4; nod++) nodindFace[nod] = incid[nod];
    gmesh->CreateGeoElement(EQuadrilateral,nodindFace,-1,index);

   }///i

  ///east
  for(int i = 0; i < n; i++){
    int incid[4] = { i*(n+1)+n, (i+1)*(n+1)+n,
                     (i+1)*(n+1)+n+(n+1)*(n+1), i*(n+1)+n+(n+1)*(n+1)
                    };
    ///face
    int index;
    TPZManVector<int,4> nodindFace(4);
    for(int nod=0; nod < 4; nod++) nodindFace[nod] = incid[nod];
    gmesh->CreateGeoElement(EQuadrilateral,nodindFace,-1,index);

   }///i

  gmesh->BuildConnectivity();

#ifdef DEBUG
{
  ofstream malhas("malhaGeo.txt");
  gmesh->Print(malhas);
}
#endif

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(3);

  TPZAutoPointer<TPZMaterial> mat = new TPZEulerEquation(1,1.4);
  cmesh->InsertMaterialObject(mat);
  TPZFMatrix val1,val2;
  cmesh->InsertMaterialObject(mat->CreateBC(mat,-1,TPZEulerEquation::EFreeSlip,val1,val2));

  cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->SetDefaultOrder(0);
  TPZCompElDisc::SetgOrder(0);

  cmesh->AutoBuild();
  TPZCompElDisc::SetTotalOrderShape(cmesh);
  TPZAutoPointer<TPZFunction> fakefunc = new TPZFakeFunction();
  for(int i = 0; i < cmesh->NElements(); i++){
    TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[i]);
    if(disc){
      disc->SetExternalShapeFunction(fakefunc);
    }
  }
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();

  return cmesh;

}

void InitialSolutionLax2D(TPZFMatrix &InitialSol, TPZCompMesh * cmesh){
  InitialSol.Redim(cmesh->NEquations(),1);
  InitialSol.Zero();
  for(int iel = 0; iel < cmesh->NElements(); iel++){
    TPZCompEl * cel = cmesh->ElementVec()[iel];
    if(!cel) continue;
    TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cel);
    if(!disc) continue;
    if(disc->NConnects() == 0) continue;
    int bl = disc->Connect(0).SequenceNumber();
    int blpos = cmesh->Block().Position(bl);
    int blocksize = cmesh->Block().Size(bl);
    
    TPZGeoEl * gel = cel->Reference();
    TPZVec<REAL> xi(3), xVec(3);
    gel->CenterPoint(gel->NSides()-1,xi);
    gel->X(xi,xVec);
    double x = xVec[0];
    double y = xVec[1];
    double rho = 1., u=0., v = 0., T = 1., w = 0.;

    if (x > 0.5 && y > 0.5)
    {
      // ZONE I
      rho = 1.;
      u = -0.75;
      v = -0.5;
    }
    
    if (x <= 0.5 && y > 0.5)
    {
      // ZONE II
      rho = 2.;
      u = -0.75;
      v = 0.5;
    }
    
    if (x <= 0.5 && y <= 0.5)
    {
      // ZONE III
      rho = 1.;
      u = 0.75;
      v = 0.5;
    }
    
    if (x > 0.5 && y <= 0.5)
    {
      // ZONE IV
      rho = 3.;
      u = 0.75;
      v = -0.5;
    }

    T=1./rho;
    const double Gamma = 1.4;
    const double Ma = 1.;
    double e = T/(Gamma*(Gamma-1.)*Ma*Ma) + 0.5*(u*u + v*v);

#ifdef ROTACAO
    ApplyRotacao(u,v,w);
#endif

    InitialSol(blpos+blocksize-20+0,0) = rho;
    InitialSol(blpos+blocksize-20+1,0) = rho*u;
    InitialSol(blpos+blocksize-20+2,0) = rho*v;
    InitialSol(blpos+blocksize-20+3,0) = rho*w;
    InitialSol(blpos+blocksize-20+4,0) = rho*e;

  }///for iel

#ifdef ROTACAO
  RotacionarMalha(cmesh);
#endif

}///method



TPZCompMesh *CreateMeshLinearConvection(int L, REAL &timeStep){
  return CreateMeshLax2D(L,timeStep);
}

/*
void InitialSolutionLinearConvection(TPZFMatrix &InitialSol, TPZCompMesh * cmesh){
  InitialSol.Redim(cmesh->NEquations(),1);
  InitialSol.Zero();
  for(int iel = 0; iel < cmesh->NElements(); iel++){
    TPZCompEl * cel = cmesh->ElementVec()[iel];
    if(!cel) continue;
    TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cel);
    if(!disc) continue;
    if(disc->NConnects() == 0) continue;
    int bl = disc->Connect(0).SequenceNumber();
    int blpos = cmesh->Block().Position(bl);
    int blocksize = cmesh->Block().Size(bl);

    TPZGeoEl * gel = cel->Reference();
    TPZVec<REAL> xi(3), xVec(3);
    gel->CenterPoint(gel->NSides()-1,xi);
    gel->X(xi,xVec);
    double x = xVec[0];
    double y = xVec[1];
    double u = 0.125;

    double xCircle = 0.25;
    double yCircle = 0.5;
    double R = 0.1;
    if( (x-xCircle)*(x-xCircle)+(y-yCircle)*(y-yCircle) <= R*R ) u = 1.;

    InitialSol(blpos+blocksize-20+0,0) = u;
    InitialSol(blpos+blocksize-20+1,0) = 0.;
    InitialSol(blpos+blocksize-20+2,0) = 0.;
    InitialSol(blpos+blocksize-20+3,0) = 0.;
    InitialSol(blpos+blocksize-20+4,0) = 0.;

  }///for iel

  TPZVec<REAL> celerity(3,0.);
  celerity[0] = 1.;
#ifdef LinearConvection
  TPZEulerEquation::SetLinearConvection(cmesh, celerity);
#endif

}//method
*/

