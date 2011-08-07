/**
 * @file
 * @brief Contains the TPZGraphElTd class which implements the graphical discontinuous triangular element.
 */
#ifndef TRIGRAPHD
#define TRIGRAPHD

#include "pzgraphel.h"
#include "pzvec.h"

/**
 * @ingroup post
 * @brief To export a graphical discontinuous triangular element. \ref post "Post processing"
 */
class TPZGraphElTd : public TPZGraphEl {
	
	public :

	TPZGraphElTd(TPZCompEl *c, TPZGraphMesh *g);
	
	virtual int NConnects();
	
	virtual int NElements();
	
	virtual MElementType Type(){return ETriangle;}
	
	virtual int ExportType(TPZDrawStyle st);
	
	virtual int NNodes();
	
	virtual TPZGraphNode *Connect(int con){ return fConnect;}
	
	virtual int NPoints(TPZGraphNode *n);
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle);
	
	virtual long EqNum(TPZVec<int> &co);
	
	protected :
	
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta);
	
	virtual void FirstIJ(int no, TPZVec<int> &co, int &incr);
	
	virtual void NextIJ(int no, TPZVec<int> &co, int incr);
	
	virtual void SetNode(int i,TPZGraphNode *gno) {
		fConnect = gno;
	}
	
	TPZGraphNode *fConnect;

};

#endif
