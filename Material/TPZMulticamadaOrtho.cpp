#include "TPZPlacaOrthotropic.h"

#include "TPZMulticamadaOrtho.h"
#include "pzmatorthotropic.h"
#include "pzgmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzbctension.h"
#include "pztempmat.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "pzskylstrmatrix.h"
// #include "pzelmat.h"
// #include "pzbndcond.h"
// #include "pzmatrix.h"
// #include "pzfmatrix.h"
// #include "pzerror"
// #include "pztempmat.h"
// #include "pzmanvector.h"
// #include <math.h>
// #include <fstream>
// using namespace std;




TPZMulticamadaOrthotropic::TPZMulticamadaOrthotropic(REAL z,REAL dx,REAL dy, int nelx, int nely){

  fGeoMesh = new TPZGeoMesh();
  fCompMesh = new TPZCompMesh(fGeoMesh);
  fZMin  =  z;
  fZMax = z;
  fNelx = nelx;
  if(! fNelx%2) fNelx++;
  fNely = nely;
  if(! fNely%2) fNely++;
  fDx = dx/fNelx;
  fDy = dy/fNely;
  fGeoMesh->NodeVec().Resize((fNelx+1)*(fNely+1));
  int ix, iy;
  TPZManVector<REAL,3> coord(3,fZMax);
  for(ix=0; ix<= nelx; ix++) {
    for(iy=0; iy<= nely; iy++) {
      coord[0] = ix*fDx-fDx*(nelx+1.)/2.;
      coord[1] = iy*fDy-fDy*(nely+1.)/2.;
      fGeoMesh->NodeVec()[ix+iy*(nelx+1)].Initialize(coord,*fGeoMesh);
    }
  }
}


void TPZMulticamadaOrthotropic::GenerateMesh(){
  fGeoMesh->BuildConnectivity();
  fCompMesh->AutoBuild();
  int nplaca = fPlacaOrth.NElements();
  int ip;
  for(ip = 0; ip<nplaca; ip++) {
    fPlacaOrth[ip].IdentifyCompEl();
  }
}


void TPZMulticamadaOrthotropic::AddPlacaOrtho(TPZMatOrthotropic *material, REAL height){

  fCompMesh->InsertMaterialObject(material);
  TPZFNMatrix<9> val1(3,3,0.),val2(3,1,0.);
  TPZBCTension *bc = new TPZBCTension(material,-material->Id(),4,val1,val2, this,fPlacaOrth.NElements());
  fCompMesh->InsertMaterialObject(bc);

  //  fPlacaOrth.Push(placa);
  int nnodes = fGeoMesh->NodeVec().NElements();
  fGeoMesh->NodeVec().Resize(nnodes+(fNelx+1)*(fNely+1));
  int ix, iy;
  TPZManVector<REAL,3> coord(3,fZMax+height);
  fZMax += height;
  for(ix=0; ix<= fNelx; ix++) {
    for(iy=0; iy<= fNely; iy++) {
      coord[0] = ix*fDx-fDx*(fNelx+1)/2;
      coord[1] = iy*fDy-fDy*(fNely+1)/2;
      fGeoMesh->NodeVec()[nnodes+ix+iy*(fNelx+1)].Initialize(coord,*fGeoMesh);
    }
  }
  int nodebase1 = nnodes - (fNelx+1)*(fNely+1);
  int elx, ely;
  TPZManVector<int,8> nodeindexes(8,-1);
  for(elx=0; elx<fNelx; elx++) {
    for(ely=0; ely<fNely; ely++) {
      nodeindexes[0] = nodebase1+elx+ely*(fNelx+1);
      nodeindexes[1] = nodebase1+elx+1+ely*(fNelx+1);
      nodeindexes[2] = nodebase1+elx+1+(ely+1)*(fNelx+1);
      nodeindexes[3] = nodebase1+elx+(ely+1)*(fNelx+1);
      int i;
      for(i=0; i<4; i++) nodeindexes[i+4] = nodeindexes[i]+(fNelx+1)*(fNely+1);
      int index;
      TPZGeoEl *gel = fGeoMesh->CreateGeoElement (ECube, nodeindexes, material->Id(), index);
      TPZGeoElBC gbc1(gel,21,-material->Id(),*fGeoMesh);
      TPZGeoElBC gbc2(gel,22,-material->Id(),*fGeoMesh);
      TPZGeoElBC gbc3(gel,23,-material->Id(),*fGeoMesh);
      TPZGeoElBC gbc4(gel,24,-material->Id(),*fGeoMesh);
      if(elx == fNelx/2 && ely == fNely/2) {
	TPZPlacaOrthotropic placa(gel,fZMax-height,fZMax);
	fPlacaOrth.Push(placa);
      }
      int nplaca = fPlacaOrth.NElements();
      if(nplaca == 0 && elx == 0 && ely == 0) {
	TPZBndCond *bc2 = new TPZBndCond(material,-100,0,val1,val2);
	fCompMesh->InsertMaterialObject(bc2);
	TPZGeoElBC gbc5(gel,0,-100,*fGeoMesh);
	val1(1,1) = 1.e12;
	val1(2,2) = 1.e12;
	TPZBndCond *bc3 = new TPZBndCond(material,-101,2,val1,val2);
	fCompMesh->InsertMaterialObject(bc3);
	TPZGeoElBC gbc6(gel,1,-101,*fGeoMesh);
	val1(1,1) = 0.;
	val1(2,2) = 1.e12;
	TPZBndCond *bc4 = new TPZBndCond(material,-102,2,val1,val2);
	fCompMesh->InsertMaterialObject(bc4);
	TPZGeoElBC gbc7(gel,2,-102,*fGeoMesh);
      }
    }
  }
}

void  TPZMulticamadaOrthotropic::Print(ostream &out){

  int i, nplaca=fPlacaOrth.NElements();
  out << "TPZMulticamadaOrthotropic::Print\n";
  out << nplaca << endl;
  for (i=0; i<nplaca; i++){
    //TPZMulticamadaOrthotropic *multcam;
    out << "placa : " << i << endl;
    fPlacaOrth[i].Print();
   //out << "quantidade de camadas :" << multcam->ZHight(placa); 
  }
}
  
REAL TPZMulticamadaOrthotropic::Height(){

  return fZMax - fZMin;
}

int TPZMulticamadaOrthotropic::NPlacas(){
  return fPlacaOrth.NElements();
}

void TPZMulticamadaOrthotropic::AnalyticTensor(TPZVec<REAL> &co, TPZFMatrix &tensor) {


  REAL z = co[2];
  REAL x = co[0];
  REAL y = co[1];
  tensor.Redim(3,3);
  REAL height = Height();
  REAL zrel = z-(fZMax-fZMin)/2.;
  if(height <= 0.) return;
  REAL height3 = height*height*height;
  tensor(0,0) += (fNX[2])/height;
  tensor(1,1) += (fNY[2])/height;
  tensor(1,0) += (fNXY[2])/height;
  tensor(0,1) += (fNXY[2])/height;
  tensor(0,0) += 12.*(fMX[2]+x*fQX[2])*zrel/(height3);
  tensor(1,1) += 12.*(fMY[2]+y*fQY[2])*zrel/height3;
  tensor(0,1) += 12.*fMXY[2]*zrel/height3;
  tensor(1,0) += 12.*fMXY[2]*zrel/height3;
  tensor(0,2) += -6.*fQX[2]*(zrel*zrel-height*height/4.)/height3;
  tensor(2,0) += -6.*fQX[2]*(zrel*zrel-height*height/4.)/height3;
  tensor(1,2) += -6.*fQY[2]*(zrel*zrel-height*height/4.)/height3;
  tensor(2,1) += -6.*fQY[2]*(zrel*zrel-height*height/4.)/height3;

}

void TPZMulticamadaOrthotropic::ComputeCenterForces() {

  int nplaca = fPlacaOrth.NElements();
  int ip;
  fMX[1] = 0.;
  fMY[1] = 0.;
  fMXY[1] = 0.;
  fNX[1] = 0.;
  fNY[1] = 0.;
  fNXY[1] = 0.;
  fQX[1] = 0.;
  fQY[1] = 0.;
  REAL zref = (fZMax+fZMin)/2.;
  TPZManVector<REAL,3> normal(3,0.),direction(3,0.);
  for(ip=0; ip<nplaca; ip++) {
    normal[0] = 1.;
    direction[0] = 1.;
    fMX[1] += fPlacaOrth[ip].Moment(zref,normal,direction);
    fNX[1] += fPlacaOrth[ip].Force(normal,direction);
    direction.Fill(0.);
    direction[1] = 1.;
    fMXY[1] += fPlacaOrth[ip].Moment(zref,normal,direction);
    fNXY[1] += fPlacaOrth[ip].Force(normal,direction);
    direction.Fill(0.);
    direction[2] = 1.;
    fQX[1] += fPlacaOrth[ip].Force(normal,direction);
    normal.Fill(0.);
    normal[1] = 1.;
    direction.Fill(0.);
    direction[1] = 1.;
    fMY[1] += fPlacaOrth[ip].Moment(zref,normal,direction);
    fNY[1] += fPlacaOrth[ip].Force(normal,direction);
    direction.Fill(0.);
    direction[2] = 1.;
    fQY[1] += fPlacaOrth[ip].Force(normal,direction);
  }
  fMX[2] = fMX[1]-fMX[0];
  fMY[2] = fMY[1]-fMY[0];
  fMXY[2] = fMXY[1]-fMXY[0];
  fNX[2] = fNX[1]-fNX[0];
  fNY[2] = fNY[1]-fNY[0];
  fNXY[2] = fNXY[1]-fNXY[0];
  fQX[2] = fQX[1]-fQX[0];
  fQY[2] = fQY[1]-fQY[0];

}

void TPZMulticamadaOrthotropic::Tensor(TPZVec<REAL> &x, int placa, TPZFMatrix &tensor) {

  REAL zmin = fPlacaOrth[placa].ZMin();
  REAL zmax = fPlacaOrth[placa].ZMax();
  REAL ksi = -1.+2.*(x[2]-zmin)/(zmax-zmin);
  TPZFNMatrix<9> tensorcomp(3,3),tensoranalytic(3,3);
  fPlacaOrth[placa].Tensor(ksi,tensorcomp);
  AnalyticTensor(x,tensoranalytic);
  tensor = tensorcomp+tensoranalytic;

}

void TPZMulticamadaOrthotropic::ComputeSolution(){
  TPZAnalysis an(fCompMesh);
  TPZSkylineStructMatrix skyl(fCompMesh);
  an.SetStructuralMatrix(skyl);
  TPZStepSolver solve;
  solve.SetDirect(ELDLt);
  an.SetSolver(solve);
  fGeoMesh->BuildConnectivity();
  fCompMesh->AutoBuild();
}
