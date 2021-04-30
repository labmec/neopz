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
	
	virtual int NConnects() override;
	
	virtual MElementType Type() override {return ETriangle;}
	
	virtual int ExportType(TPZDrawStyle st) override;
	
	virtual int NNodes() override;
	
	virtual int NElements() override;
	
	virtual int NPoints(TPZGraphNode *n) override;
	
	virtual TPZGraphNode *Connect(int64_t i) override {return fConnects[i];}
	
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
		fConnects[i] = gno;
	}
	
	/** @brief Graphical nodes vector (by connect of the computational element) */
	TPZGraphNode *fConnects[7];
	
};

#endif
