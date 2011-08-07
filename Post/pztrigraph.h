/**
 * @file
 * @brief Contains the TPZGraphElT class which implements the graphical triangular element.
 */
#ifndef TRIGRAPH
#define TRIGRAPH

#include "pzgraphel.h"
#include "pzvec.h"

/**
 * @ingroup post
 * @brief To export a graphical triangular element. \ref post "Post processing"
 */
class TPZGraphElT : public TPZGraphEl {
	
	public :
	
	
	TPZGraphElT(TPZCompEl *c, TPZGraphMesh *g);
	
	virtual int NConnects();
	
	virtual MElementType Type() {return ETriangle;}
	
	virtual int ExportType(TPZDrawStyle st);
	
	virtual int NNodes();
	
	virtual int NElements();
	
	virtual int NPoints(TPZGraphNode *n);
	
	virtual TPZGraphNode *Connect(int i) {return fConnects[i];}
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle);
	
	virtual long EqNum(TPZVec<int> &co);
	
	protected :
	
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta);
	
	virtual void FirstIJ(int no, TPZVec<int> &co, int &incr);
	
	virtual void NextIJ(int no, TPZVec<int> &co, int incr);
	
	virtual void SetNode(int i,TPZGraphNode *gno) {
		fConnects[i] = gno;
	}
	
	TPZGraphNode *fConnects[7];

};

#endif
