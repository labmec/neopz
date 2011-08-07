// TPZGraphElQ3dd.h: interface for the TPZGraphElQ3dd class.
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
	
	TPZGraphElQ3dd(TPZCompEl *cel, TPZGraphMesh *gmesh);
	
	virtual ~TPZGraphElQ3dd(void);
	
	virtual int NConnects(){ return 1;}
	
	virtual MElementType Type(){return ECube;}
	
	virtual int ExportType(TPZDrawStyle st);
	
	virtual int NNodes();
	
	virtual TPZGraphNode *Connect(int i);
	
	virtual int NPoints(TPZGraphNode *n);
	
	virtual int NElements();
	
	virtual	void SetNode(int i,TPZGraphNode *n);
	
	virtual long EqNum(TPZVec<int> &co);
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle);
	
	
protected:
	
	virtual void FirstIJ(int connect, TPZVec<int> &co, int &incr);
	
	virtual void NextIJ(int connect, TPZVec<int> &co, int incr);
	
	protected :   
	
	TPZGraphNode *fConnect;

};

#endif // !defined(AFX_TPZGRAPHELQ3DD_H__4DDADE46_92E7_11D4_B7FB_00500464279E__INCLUDED_)
