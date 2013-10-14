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
	
	
	/** @brief Constructor for graphical element to computational triangular element */
	TPZGraphElT(TPZCompEl *c, TPZGraphMesh *g);
	
	virtual int NConnects();
	
	virtual MElementType Type() {return ETriangle;}
	
	virtual int ExportType(TPZDrawStyle st);
	
	virtual int NNodes();
	
	virtual int NElements();
	
	virtual int NPoints(TPZGraphNode *n);
	
	virtual TPZGraphNode *Connect(long i) {return fConnects[i];}
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle);
	
	virtual long EqNum(TPZVec<int> &co);
	
	protected :
	
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta);
	
	virtual void FirstIJ(int no, TPZVec<int> &co, int &incr);
	
	virtual void NextIJ(int no, TPZVec<int> &co, int incr);
	
	virtual void SetNode(long i,TPZGraphNode *gno) {
		fConnects[i] = gno;
	}
	
	/** @brief Graphical nodes vector (by connect of the computational element) */
	TPZGraphNode *fConnects[7];
	
};

#endif
