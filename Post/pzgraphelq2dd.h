/**
 * @file
 * @brief Contains the TPZGraphElQ2dd class which implements the graphical two-dimensional discontinuous element.
 */

#ifndef PZGRAPHELQ2D
#define PZGRAPHELQ2D

#include "pzgraphel.h"
#include "pzvec.h"

/**
 * @ingroup post
 * @brief To export a graphical two-dimensional discontinuous element. \ref post "Post processing"
 */
class TPZGraphElQ2dd : public TPZGraphEl {
public:
	/** @brief Constructor for graphical element to computational quadrilateral discontinuous element */
	TPZGraphElQ2dd(TPZCompEl *cel, TPZGraphMesh *gmesh);
	
	virtual ~TPZGraphElQ2dd(void);
	
	virtual int NConnects() override { return 1;}
	
	virtual MElementType Type() override {return EQuadrilateral;}
	
	virtual int ExportType(TPZDrawStyle st) override;
	
	virtual int NNodes() override;
	
	virtual TPZGraphNode *Connect(int64_t i) override;
	
	virtual int NPoints(TPZGraphNode *n) override;
	
	virtual int NElements() override;
	
	virtual	void SetNode(int64_t i,TPZGraphNode *n) override;
	
	virtual int64_t EqNum(TPZVec<int> &co) override;
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle) override;
    
    /** @brief the parametric dimension of the element */
    virtual int Dimension() override
    {
        return 2;
    }

	
	
protected:
	
	virtual void FirstIJ(int connect, TPZVec<int> &co, int &incr) override;
	
	virtual void NextIJ(int connect, TPZVec<int> &co, int incr) override;
	
	protected :   
	
	/** @brief Graphical node (connect) to discontinuous graphical element */
	TPZGraphNode *fConnect;
	
};

#endif
