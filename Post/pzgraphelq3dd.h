/**
 * @file
 * @brief Contains the TPZGraphElQ3dd class which implements the graphical three dimensional discontinuous element.
 */

#if !defined(AFX_TPZGRAPHELQ3DD_H__4DDADE46_92E7_11D4_B7FB_00500464279E__INCLUDED_)
#define AFX_TPZGRAPHELQ3DD_H__4DDADE46_92E7_11D4_B7FB_00500464279E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "pzgraphel.h"
#include "pzvec.h"

/**
 * @ingroup post
 * @brief To export a graphical three dimensional discontinuous element. \ref post "Post processing"
 */
class TPZGraphElQ3dd : public TPZGraphEl{
public:
	
	/** @brief Constructor for graphical element to computational hexahedra discontinuous element */
	TPZGraphElQ3dd(TPZCompEl *cel, TPZGraphMesh *gmesh);
	
	virtual ~TPZGraphElQ3dd(void);
	
	virtual int NConnects() override { return 1;}
	
	virtual MElementType Type() override {return ECube;}
	
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
        return 3;
    }

	
	
protected:
	
	virtual void FirstIJ(int connect, TPZVec<int> &co, int &incr) override;
	
	virtual void NextIJ(int connect, TPZVec<int> &co, int incr) override;
	
	protected :   
	
	/** @brief Graphical node (connect) to discontinuous graphical element */
	TPZGraphNode *fConnect;
	
};

#endif // !defined(AFX_TPZGRAPHELQ3DD_H__4DDADE46_92E7_11D4_B7FB_00500464279E__INCLUDED_)
