
#ifndef ITERATIVEANALYSISH
#define ITERATIVEANALYSISH

#include <iostream>
#include "pzanalysis.h"
#include "pzreal.h"
using namespace std;
class TPZCompMesh;
class TPZCompEl;
class TPZMaterial;
/******************/

class TPZIterativeAnalysis : public TPZAnalysis {

  clock_t fBegin,fInit;

public:

  TPZIterativeAnalysis(TPZCompMesh *mesh,std::ostream &out = cout);

  void IterativeProcess(ostream &out,REAL tol,int numiter,TPZMaterial *mat,int marcha=1,int resolution=0);

  void IterativeProcessTest(ostream &out,REAL tol,int numiter,TPZMaterial *mat,int marcha,int resolution=0);

  void CoutTime(clock_t &start,char *title);

  virtual void SetDeltaTime(TPZCompMesh *CompMesh,TPZMaterial *mat);
};

#endif

