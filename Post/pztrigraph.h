#ifndef TRIGRAPH
#define TRIGRAPH

#include "pzgraphel.h"
#include "pzvec.h"

class TPZGraphElT : public TPZGraphEl {

public :


TPZGraphElT(TPZCompEl *c, TPZGraphMesh *g);

virtual int NConnects();

 virtual MElementType Type() {return ETriangle;}

virtual int NElements();

virtual int NPoints(TPZGraphNode *n);

virtual TPZGraphNode *Connect(int i) {return fConnects[i];}

virtual void Connectivity(TPZDrawStyle st = EDXStyle);

virtual long EqNum(TPZVec<int> &co);

protected :

REAL QsiEta(int i, int imax);

virtual void FirstIJ(int no, TPZVec<int> &co, int &incr);

virtual void NextIJ(int no, TPZVec<int> &co, int incr);

virtual void SetNode(int i,TPZGraphNode *gno) {
	fConnects[i] = gno;
}

	TPZGraphNode *fConnects[7];


};

#endif
