/**
 * @file
 * @brief Contains the TPZGraphElQ2d class which implements the graphical two dimensional element.
 */
#ifndef PZGRAPHELQ2D
#define PZGRAPHELQ2D

#include "pzgraphel.h"
#include "pzvec.h"

/**
 * @ingroup post
 * @brief To export a graphical two dimensional element. \ref post "Post processing"
 */
class TPZGraphElQ2d : public TPZGraphEl {
public:
	/** @brief Constructor for graphical element to computational quadrilateral element */
	TPZGraphElQ2d(TPZCompEl *cel, TPZGraphMesh *gmesh);
	
	virtual ~TPZGraphElQ2d(void);
	
	virtual int NConnects() override { return 9;}
	
	virtual MElementType Type() override { return EQuadrilateral;}
	
	virtual int ExportType(TPZDrawStyle st) override;
	
	virtual int NNodes() override;
	
	virtual TPZGraphNode *Connect(int64_t i) override;
	
	virtual int NPoints(TPZGraphNode *n) override;
	
	virtual int NElements() override;
	
	virtual void SetNode(int64_t i,TPZGraphNode *n) override;
	
	virtual int64_t EqNum(TPZVec<int> &co) override;
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle) override;

    /** @brief the parametric dimension of the element */
    virtual int Dimension() override
    {
        return 2;
    }

	
protected:
	
	virtual void FirstIJ(int connect, TPZVec<int> &co, int &incr) override;
	
	virtual void NextIJ(int connect,TPZVec<int> &co, int incr) override;
	
	protected :   
	
	/** @brief Graphical nodes vector (by connect of the computational element) */
	TPZGraphNode *fConnects[9];	
	
};

#endif
