#ifndef GRAFEL1DDH
#define GRAFEL1DDH


#include "pzgraphel.h"
#include "pzvec.h"

class TPZGraphMesh;
class TPZGraphNode;
class TPZBlock;

class TPZGraphEl1dd : public TPZGraphEl
{
public:
	TPZGraphEl1dd(TPZCompEl *ce, TPZGraphMesh *gg);

virtual int NNodes();

virtual int NPoints(TPZGraphNode *n);

virtual int NElements();

virtual void Connectivity(TPZDrawStyle st = EDXStyle);

void Print(ostream &out);

virtual long EqNum(TPZVec<int> &co);


protected:

virtual void FirstIJ(int no, TPZVec<int> &co, int &incr);

virtual void NextIJ(int no, TPZVec<int> &co, int incr);

//	void ComputeSequence(TPZGraphNode *n,int *ibound,int *incr);

	TPZGraphNode *fConnect;

	virtual void SetNode(int i,TPZGraphNode *gno) {
		fConnect = gno;
	}

/** Jorge 8/6/2001 */
	int NConnects() { return 1; }
	MElementType Type() {return EOned;}
	TPZGraphNode *Connect(int i) {return fConnect;}

};

#endif
