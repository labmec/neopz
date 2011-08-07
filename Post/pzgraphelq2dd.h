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
	
	TPZGraphElQ2dd(TPZCompEl *cel, TPZGraphMesh *gmesh);
	
	virtual ~TPZGraphElQ2dd(void);
	
	virtual int NConnects(){ return 1;}
	
	virtual MElementType Type(){return EQuadrilateral;}
	
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

#endif
