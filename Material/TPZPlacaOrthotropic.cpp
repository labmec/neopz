#include "pzintel.h"
#include "TPZPlacaOrthotropic.h"
#include "pzgnode.h"
#include "pzgeoel.h"
#include "pzquad.h"

TPZPlacaOrthotropic::TPZPlacaOrthotropic(TPZGeoEl *gel,REAL zmin, REAL zmax ){ 

  fGeoEl = gel;
  fH = zmax-zmin;
  fZMin = zmin;
  fZMax = zmax;
  fIntel = dynamic_cast<TPZInterpolatedElement *>(gel->Reference());
  if(fIntel) {
    fTensorVar = fIntel->Material()->VariableIndex("Tensor");
  }
}

TPZPlacaOrthotropic::TPZPlacaOrthotropic(){ 

  fGeoEl = 0;
  fH = 0.;
  fZMin = 0.;
  fZMax = 0.;
  fIntel = 0;
  fTensorVar = -1;
}

void TPZPlacaOrthotropic::Tensor(REAL ksi, TPZFMatrix &T) {
  if(fTensorVar == -1) {
    if(!fIntel) fIntel = dynamic_cast<TPZInterpolatedElement *>(fGeoEl->Reference());
    if(!fIntel) return;
    fTensorVar = fIntel->Material()->VariableIndex("Tensor");
  }
  if(fTensorVar == -1) return;
  //  REAL ksi = -1+2.*(z-fZMin)/(fZMax-fZMin);
  if(ksi < -1. || ksi > 1.) return;
  TPZManVector<REAL,3> co(3,0.);
  co[2] = ksi;
  TPZManVector<REAL> tensor(9);
  fIntel->Solution(co,fTensorVar,tensor);
  int i;
  for(i=0; i<9; i++) {
    T(i%3,i/3) = tensor[i];
  }
}

REAL TPZPlacaOrthotropic::Moment(REAL zref, TPZVec<REAL> &normal, TPZVec<REAL> &direction){
  TPZInt1d rule(4);
  int np = rule.NPoints();
  TPZManVector<REAL,3> pos(1,0.);
  TPZFNMatrix<9> tensor(3,3,0.);
  int ip;
  REAL moment = 0.;
  for(ip = 0; ip<np; ip++) {
    REAL weight;
    rule.Point(ip,pos,weight);
    Tensor(pos[0],tensor);
    REAL z = fZMin + (fZMax-fZMin)*(pos[0]+1.)/2.;
    REAL tension = 0.;
    int n,d;
    for(n=0; n<3; n++) {
      for(d=0; d<3; d++) {
	tension += normal[n]*tensor(n,d)*direction[d];
      }
    }
    moment += weight*tension*(z-zref);
  }
  moment *= fH/2.;
  return moment;

}

REAL TPZPlacaOrthotropic::Force(TPZVec<REAL> &normal, TPZVec<REAL> &direction){
  TPZInt1d rule(4);
  int np = rule.NPoints();
  TPZManVector<REAL,3> pos(1,0.);
  TPZFNMatrix<9> tensor(3,3,0.);
  int ip;
  REAL force = 0.;
  for(ip = 0; ip<np; ip++) {
    REAL weight;
    rule.Point(ip,pos,weight);
    Tensor(pos[0],tensor);
    REAL tension = 0.;
    int n,d;
    for(n=0; n<3; n++) {
      for(d=0; d<3; d++) {
	tension += normal[n]*tensor(n,d)*direction[d];
      }
    }
    force += weight*tension;
  }
  force *= fH/2.;
  return force;
}


// TPZPlacaOrthotropic::CriaNos(int num, TPZGeomesh &geomesh, double list[20][3]){
//   geomesh.NodeVec().Resize(num);
//   TPZVec<REAL> coord(3);
//   int i;
//   for (i=0; i<num; i++){
//     coord[0] = list[i][0];
//     coord[1] = list[i][1];
//     coord[2] = list[i][2];
//     geomesh.NodeVec()[i].Initialize(coord, geomesh);
//   }
// }

//  orto = new TPZMatOrthotropic(1, naxes, 69.e06, 69.e06, 69.e06, 0.33, 0.33, 0.33, 5.15378e06, 5.15378e06, 5.15378e06) ;
//  TPZMaterial *orto2 = new TPZMatOrthotropic(1, naxes, 69.e06, 69.e06, 69.e06, 0.33, 0.33, 0.33, 5.15378e06, 5.15378e06, 5.15378e06) ;

void TPZPlacaOrthotropic::Print(){

  fIntel->Print();
  cout << "Plate height : " ;
  cout << this->Height();//this é o objeto placa

}

//se fosse TPZMulticamadaOrthotropic::Print, o this seria multcam
//se a funcão fosse declarada como void Print(){...}, não haveria This.
//caso chamemos a funcão Print com o objeto placa2, ao invés de placa->Print, o This é o objeto placa2.
  
  
