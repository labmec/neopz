#ifndef GRAFEL1DH
#define GRAFEL1DH


#include "pzgraphel.h"
#include "pzvec.h"

class TPZGraphMesh;
class TPZGraphNode;
class TPZBlock;

class TPZGraphEl1d : public TPZGraphEl
{
public:
	TPZGraphEl1d(TPZCompEl *ce, TPZGraphMesh *gg);

virtual int NConnects() { return 3;}

 virtual MElementType Type() {return EOned;}

virtual TPZGraphNode *Connect(int i) {return fConnects[i];}

virtual int NNodes();

virtual int NPoints(TPZGraphNode *n);

virtual int NElements();

virtual void Connectivity(TPZDrawStyle st = EDXStyle);

void Print(ostream &out);

virtual long EqNum( TPZVec<int> &co);


protected:

virtual void FirstIJ(int no, TPZVec<int> &co, int &incr);

virtual void NextIJ(int no,  TPZVec<int> &co, int incr);

//	void ComputeSequence(TPZGraphNode *n,int *ibound,int *incr);

	TPZGraphNode *fConnects[3];

virtual void SetNode(int i,TPZGraphNode *gno) {
	fConnects[i] = gno;
}


};

#endif
