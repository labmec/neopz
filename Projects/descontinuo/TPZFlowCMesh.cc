
#include "pzcompel.h"
#include "pzgeoel.h"
#include "TPZConservationLaw.h"
#include "TPZFlowCMesh.h"
#include "TPZAgglomerateEl.h"

TPZFlowCompMesh::~TPZFlowCompMesh(){

}

void TPZFlowCompMesh::NotCreate(){
  fExists = 0;
}

int TPZFlowCompMesh::Exists(){
  return fExists;
}

TPZFlowCompMesh::TPZFlowCompMesh(TPZGeoMesh* gr) : TPZCompMesh(gr) {
  fExists = 1;
}


REAL TPZFlowCompMesh::MaxVelocityOfMesh(int nstate,REAL gamma){

  int nel = ElementVec().NElements(),i;
  TPZManVector<REAL> density(1),sol(nstate),velocity(1);
  REAL maxvel = 0.0,veloc,sound,press;
  TPZVec<REAL> param(3,0.);
  for(i=0;i<nel;i++){
    TPZCompEl *com = ElementVec()[i];
    if(!com) continue;
    int type = com->Type();
    if(type == EInterface) continue;
    TPZMaterial *mat = com->Material();
    if(!mat){
      cout << "TPZFlowCompMesh::MaxVelocityOfMesh ERROR: null material\n";
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
      cout << "TPZFlowCompMesh::MaxVelocityOfMesh unknown case\n";
    }
    com->Solution(param,1,density);
    if(density[0] < 0.0){
      cout << "TPZFlowCompMesh::MaxVelocityOfMesh minus density\n";
    }
    com->Solution(param,6,velocity);
    com->Solution(param,5,sol);
    TPZConservationLaw *law = dynamic_cast<TPZConservationLaw *>(mat);
    press = law->Pressure(sol);
    if(press < 0.0){
      cout << "TPZFlowCompMesh::MaxVelocityOfMesh minus pressure\n";
    }
    sound = sqrt(gamma*press/density[0]);
    veloc = velocity[0] + sound;
    if(veloc > maxvel) maxvel = veloc;
  }
  return maxvel;

}


