
#ifndef ITERATIVEANALYSISH
#define ITERATIVEANALYSISH

#include <iostream>
#include <stdio.h>
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

  TPZIterativeAnalysis(TPZCompMesh *mesh,ostream &out = cout);

	void IterativeProcess(std::string &name,REAL tol,int numiter,TPZAutoPointer<TPZMaterial> mat,int marcha=1,int resolution=0);

  void IterativeProcessTest(std::string &name,REAL tol,int numiter,TPZAutoPointer<TPZMaterial> mat,int marcha,int resolution=0);

  void CoutTime(clock_t &start,const char *title);

  virtual void SetDeltaTime(TPZCompMesh *CompMesh,TPZMaterial *mat);

  void ResetReference(TPZCompMesh *aggcmesh);

  void SetReference(TPZCompMesh *aggcmesh);
};

#endif

