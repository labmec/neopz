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
	
	/** @brief Constructor for graphical element to computational triangular discontinuous element */
	TPZGraphElTd(TPZCompEl *c, TPZGraphMesh *g);
	
	virtual int NConnects();
	
	virtual int NElements();
	
	virtual MElementType Type(){return ETriangle;}
	
	virtual int ExportType(TPZDrawStyle st);
	
	virtual int NNodes();
	
	virtual TPZGraphNode *Connect(int64_t con){ return fConnect;}
	
	virtual int NPoints(TPZGraphNode *n);
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle);
	
	virtual int64_t EqNum(TPZVec<int> &co);
    
    /** @brief the parametric dimension of the element */
    virtual int Dimension()
    {
        return 2;
    }

	
	protected :
	
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta);
	
	virtual void FirstIJ(int no, TPZVec<int> &co, int &incr);
	
	virtual void NextIJ(int no, TPZVec<int> &co, int incr);
	
	virtual void SetNode(int64_t i,TPZGraphNode *gno) {
		fConnect = gno;
	}
	
	/** @brief Graphical node (connect) to discontinuous graphical element */
	TPZGraphNode *fConnect;
	
};

#endif
