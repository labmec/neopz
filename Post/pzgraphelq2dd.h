#ifndef PZGRAPHELQ2D
#define PZGRAPHELQ2D

#include "pzgraphel.h"
#include "pzvec.h"

class TPZGraphElQ2dd : public TPZGraphEl {
public:

TPZGraphElQ2dd(TPZCompEl *cel, TPZGraphMesh *gmesh);

virtual ~TPZGraphElQ2dd(void);

virtual int NConnects(){ return 1;}

 virtual MElementType Type(){return EQuadrilateral;}

virtual TPZGraphNode *Connect(int i);

virtual int NPoints(TPZGraphNode *n);

virtual int NElements();

virtual	void SetNode(int i,TPZGraphNode *n);

virtual long EqNum(TPZVec<int> &co);

virtual void Connectivity(TPZDrawStyle st = EDXStyle);


protected:

virtual void FirstIJ(int connect, TPZVec<int> &co, int &incr);

virtual void NextIJ(int connect, TPZVec<int> &co, int incr);

//	void ComputeSequence(TPZGraphNode *n,int *ibound,int *incr);

protected :   

	TPZGraphNode *fConnect;


};

#endif
