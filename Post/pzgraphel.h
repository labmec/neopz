#ifndef GRAFELH
#define GRAFELH


#include <iostream>
using namespace std;
#include <string.h>
#include "pzcompel.h"
#include "pzgraphnode.h"
#include "pzgraphmesh.h"
#include "pzvec.h"

class TPZGraphMesh;
class TPZGraphNode;
class TPZBlock;

class TPZGraphEl
{
public:

  TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode **connectvec);

  TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode *&connect);

  virtual ~TPZGraphEl(void);

  virtual int NConnects() = 0;

  int Id() {return fId;}

  virtual MElementType Type() = 0;

  virtual TPZGraphNode *Connect(int con) = 0;

  void SetId(int id) { fId = id;}

  virtual int NPoints(TPZGraphNode *n) = 0;

  virtual int NElements() = 0;

  virtual void SetNode(int i,TPZGraphNode *n);

  virtual void Connectivity(TPZDrawStyle st = EDXStyle) = 0;

  void DrawCo(TPZGraphNode *n, TPZDrawStyle st);

  void DrawSolution(TPZGraphNode *n,TPZBlock &Sol, TPZDrawStyle st);
  void DrawSolution(TPZGraphNode *n,int solind, TPZDrawStyle st);
  void DrawSolution(TPZGraphNode *n,TPZVec<int> &solind, TPZDrawStyle st);

  void Print(ostream &out);

  virtual long EqNum(TPZVec<int> &co) = 0;


protected:
	TPZCompEl *fCompEl;
	TPZGraphMesh *fGraphMesh;

virtual void FirstIJ(int connect, TPZVec<int> &co, int &incr) = 0;

virtual void NextIJ(int connect, TPZVec<int> &co, int incr) = 0;

//virtual void ComputeSequence(TPZGraphNode *n,int *ibound,int *incr) = 0;
	int ConnectNum(TPZGraphNode *n);

protected :   
	int fId;

virtual REAL QsiEta(int i, int imax);

};

#endif
