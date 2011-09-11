
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzconslaw.h"
#include "TPZFlowCMesh.h"
#include "TPZAgglomerateEl.h"

#include <stdio.h>
#include <iostream>
using namespace std;

TPZFlowCompMesh1::~TPZFlowCompMesh1() {
}

void TPZFlowCompMesh1::NotCreate() {
  fExists = 0;
}

int TPZFlowCompMesh1::Exists() {
  return fExists;
}

TPZFlowCompMesh1::TPZFlowCompMesh1(TPZGeoMesh* gr) : TPZCompMesh(gr) {
  fExists = 1;
}


REAL TPZFlowCompMesh1::MaxVelocityOfMesh(int nstate,REAL gamma) {

  int nel = ElementVec().NElements(),i;
  TPZManVector<REAL> density(1),sol(nstate),velocity(1);
  REAL maxvel = 0.0,veloc,sound,press;
  TPZVec<REAL> param(3,0.);
  for(i=0;i<nel;i++){
    TPZCompEl *com = ElementVec()[i];
    if(!com) continue;
    int type = com->Type();
    if(type == EInterface) continue;
    TPZMaterial* mat = com->Material().operator->();
    if(!mat){
      cout << "TPZFlowCompMesh1::MaxVelocityOfMesh ERROR: null material\n";
      continue;
    }
    if(mat->Id() < 0) continue;//BC element
    TPZGeoEl *geo = com->Reference();
    if(type == EAgglomerate || type == EDiscontinuous){
      //ponto deformado
      dynamic_cast<TPZCompElDisc *>(com)->InternalPoint(param);
    } else if(geo){
      //interpolados: ponto no elemento mestre
      geo->CenterPoint(geo->NSides()-1,param);//com->Solution(sol,j+100,sol2);
    } else {
      cout << "TPZFlowCompMesh1::MaxVelocityOfMesh unknown case\n";
    }
    com->Solution(param,1,density);
    if(density[0] < 0.0){
      cout << "TPZFlowCompMesh1::MaxVelocityOfMesh minus density\n";
    }
    com->Solution(param,6,velocity);
    com->Solution(param,5,sol);
   TPZConservationLaw *law = dynamic_cast< TPZConservationLaw *> (mat);
   // press = law->Pressure(sol);
	//  TPZAutoPointer<TPZConservationLaw> law = (dynamic_cast<TPZAutoPointer<TPZConservationLaw> >(mat));
	  press = law->Pressure(sol);
    if(press < 0.0){
      cout << "TPZFlowCompMesh1::MaxVelocityOfMesh minus pressure\n";
    }
    sound = sqrt(gamma*press/density[0]);
    veloc = velocity[0] + sound;
    if(veloc > maxvel) maxvel = veloc;
  }
  return maxvel;

}


