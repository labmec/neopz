#include "TPZPlacaOrthotropic.h"

#include "TPZMulticamadaOrtho.h"
#include "pzmatorthotropic.h"
#include "pzgmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pzbctension.h"
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
  fDx = dx;
  fDy = dy;
  fNelx = nelx;
  fNely = nely;
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
  fCompMesh->AutoBuild();
}


void TPZMulticamadaOrthotropic::AddPlacaOrtho(TPZMatOrthotropic *material, REAL height){

  fCompMesh->InsertMaterialObject(material);
  TPZFNMatrix<9> val1(3,3,0.),val2(3,1,0.);
  TPZBCTension *bc = new TPZBCTension(material,-material->Id(),4,val1,val2, this);
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
	TPZPlacaOrthotropic placa(gel,height);
	fPlacaOrth.Push(placa);
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



  

