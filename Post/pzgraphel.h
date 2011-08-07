/**
 * @file
 * @brief Contains the TPZGraphEl class which implements the graphical one-, two- and three-dimensional element.
 */
#ifndef GRAFELH
#define GRAFELH

#include "pzcompel.h"
#include "pzgraphnode.h"
#include "pzgraphmesh.h"
#include "pzvec.h"

#include <iostream>
class TPZGraphMesh;
class TPZGraphNode;
class TPZBlock;

/**
 * @ingroup post
 * @brief Abstract class to graphical one-, two- and three-dimensional element. \ref post "Post processing"
 */
class TPZGraphEl
{
public:
	
	TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode **connectvec);
	
	TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode *&connect);
	
	virtual ~TPZGraphEl(void);
	
	virtual int NConnects() = 0;
	
	int Id() {return fId;}
	
	virtual MElementType Type() = 0;
	
	virtual int ExportType(TPZDrawStyle st) = 0;
	
	virtual int NNodes() = 0;
	
	virtual TPZGraphNode *Connect(int con) = 0;
	
	void SetId(int id) { fId = id;}
	
	virtual int NPoints(TPZGraphNode *n) = 0;
	
	virtual int NElements() = 0;
	
	virtual void SetNode(int i,TPZGraphNode *n);
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle) = 0;
	
	void DrawCo(TPZGraphNode *n, TPZDrawStyle st);
	
	void DrawSolution(TPZGraphNode *n,TPZBlock &Sol, TPZDrawStyle st);
	void DrawSolution(TPZGraphNode *n,int solind, TPZDrawStyle st);
	void DrawSolution(TPZGraphNode *n,TPZVec<int> &solind, TPZDrawStyle st);
	
	void Print(std::ostream &out);
	
	virtual long EqNum(TPZVec<int> &co) = 0;
	
	
protected:
	TPZCompEl *fCompEl;
	TPZGraphMesh *fGraphMesh;
	
	virtual void FirstIJ(int connect, TPZVec<int> &co, int &incr) = 0;
	
	virtual void NextIJ(int connect, TPZVec<int> &co, int incr) = 0;
	
	int ConnectNum(TPZGraphNode *n);
	
	protected :   
	int fId;
	/**
	 * @brief This method maps the index of a point to parameter space as a function
	 * of the number of divisions
	 */
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta);
	
};

#endif
