
#ifndef TPZFLOWCOMPMESHPP
#define TPZFLOWCOMPMESHPP

class TPZMaterial;
class TPZGeoMesh;
#include <iostream>
#include "pzcmesh.h"
/*******************/
class TPZFlowCompMesh1 : public TPZCompMesh {

  int fExists;

 public:

  //TPZFlowCompMesh1();/**empty constructor*/

  TPZFlowCompMesh1(TPZGeoMesh* gr);
  
  ~TPZFlowCompMesh1();

  REAL MaxVelocityOfMesh(int nstate,REAL gamma);

  void SetDeltaTime(TPZMaterial *mat);

  int Exists();

  void NotCreate();
};

#endif

