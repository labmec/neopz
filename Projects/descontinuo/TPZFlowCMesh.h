#include "pzcompel.h"
#include "pzgeoel.h"
#include "TPZConservationLaw.h"


#ifndef TPZFLOWCOMPMESHPP
#define TPZFLOWCOMPMESHPP

class TPZMaterial;
class TPZGeoMesh;
#include <iostream>
#include "pzcmesh.h"
using namespace std;
/*******************/
class TPZFlowCompMesh : public TPZCompMesh {

 public:

  //TPZFlowCompMesh();/**empty constructor*/

  TPZFlowCompMesh(TPZGeoMesh* gr);
  
  ~TPZFlowCompMesh(){};

  REAL MaxVelocityOfMesh(int nstate,REAL gamma);

  void SetDeltaTime(TPZMaterial *mat);
};



inline TPZFlowCompMesh::TPZFlowCompMesh(TPZGeoMesh* gr) : TPZCompMesh(gr) {

}


inline REAL TPZFlowCompMesh::MaxVelocityOfMesh(int nstate,REAL gamma){

  int nel = ElementVec().NElements(),i;
  TPZManVector<REAL> density(1),sol(nstate),velocity(1);
  REAL maxvel = 0.0,veloc,sound,press;
  TPZVec<REAL> param(3,0.);
  for(i=0;i<nel;i++){
    TPZCompEl *com = ElementVec()[i];
    if(!com) continue;
    TPZMaterial *mat = com->Material();
    TPZGeoEl *geo = com->Reference();
    if(!mat || !geo){
      cout << "TPZFlowCompMesh::MaxVelocityOfMesh ERROR: null material or element\n";
      continue;
    }
    int dim = mat->Dimension();
    int dimel = geo->Dimension();
    if(dimel != dim) continue;
    geo->CenterPoint(geo->NSides()-1,param);//com->Solution(sol,j+100,sol2);
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

#endif

