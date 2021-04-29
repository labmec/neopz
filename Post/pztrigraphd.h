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
	
	virtual int NConnects() override;
	
	virtual int NElements() override;
	
	virtual MElementType Type() override {return ETriangle;}
	
	virtual int ExportType(TPZDrawStyle st) override;
	
	virtual int NNodes() override;
	
	virtual TPZGraphNode *Connect(int64_t con) override{ return fConnect;}
	
	virtual int NPoints(TPZGraphNode *n) override;
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle) override;
	
	virtual int64_t EqNum(TPZVec<int> &co) override;
    
    /** @brief the parametric dimension of the element */
    virtual int Dimension() override
    {
        return 2;
    }

	
	protected :
	
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta) override;
	
	virtual void FirstIJ(int no, TPZVec<int> &co, int &incr) override;
	
	virtual void NextIJ(int no, TPZVec<int> &co, int incr) override;
	
	virtual void SetNode(int64_t i,TPZGraphNode *gno) override {
		fConnect = gno;
	}
	
	/** @brief Graphical node (connect) to discontinuous graphical element */
	TPZGraphNode *fConnect;
	
};

#endif
