
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

  TPZFlowCompMesh(TPZGeoMesh* gr);
  
  ~TPZFlowCompMesh(){};

  REAL MaxVelocityOfMesh(int nstate,REAL gamma);

  void SetDeltaTime(TPZMaterial *mat);
};

#endif

