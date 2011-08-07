/**
 * @file
 * @brief Contains the TPZGraphEl1d class which implements the graphical one dimensional element.
 */
#ifndef GRAFEL1DH
#define GRAFEL1DH

#include "pzgraphel.h"
#include "pzvec.h"

class TPZGraphMesh;
class TPZGraphNode;
class TPZBlock;

/**
 * @ingroup post
 * @brief To export a graphical one dimensional element. \ref post "Post processing"
 */
class TPZGraphEl1d : public TPZGraphEl
{
public:
	TPZGraphEl1d(TPZCompEl *ce, TPZGraphMesh *gg);
	
	virtual int NConnects() { return 3;}
	
	virtual MElementType Type() {return EOned;}
	
	virtual int ExportType(TPZDrawStyle st);
	
	virtual int NNodes();
	
	virtual TPZGraphNode *Connect(int i) {return fConnects[i];}
	
	virtual int NPoints(TPZGraphNode *n);
	
	virtual int NElements();
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle);
	
	void Print(std::ostream &out);
	
	virtual long EqNum( TPZVec<int> &co);
	
	
protected:
	
	virtual void FirstIJ(int no, TPZVec<int> &co, int &incr);
	
	virtual void NextIJ(int no,  TPZVec<int> &co, int incr);
	
	TPZGraphNode *fConnects[3];
	
	virtual void SetNode(int i,TPZGraphNode *gno) {
		fConnects[i] = gno;
	}
	
};

#endif
