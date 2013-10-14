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
template <class>
class TPZBlock;

/**
 * @ingroup post
 * @brief Abstract class to graphical one-, two- and three-dimensional element. \ref post "Post processing"
 */
class TPZGraphEl
{
public:
	/** @brief Constructor of the graphical element */
	TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode **connectvec);
	/** @brief Constructor of the graphical element */	
	TPZGraphEl(TPZCompEl *cel, TPZGraphMesh *gmesh, TPZGraphNode *&connect);
	/** @brief Default destructor */
	virtual ~TPZGraphEl(void);
	/** @brief Number of connects for the element */
	virtual int NConnects() = 0;
	/** @brief Get the Id of the graphical element */
	long Id() {return fId;}
	/** @brief Get the type of the graphical element */
	virtual MElementType Type() = 0;
	
	/** @brief Sets the style to export (format) */
	virtual int ExportType(TPZDrawStyle st) = 0;
	/** @brief Number of corner nodes (geometric information) */
	virtual int NNodes() = 0;
	
	/** @brief Return the graphical connect */
	virtual TPZGraphNode *Connect(long con) = 0;
	/** @brief Set graphical element id */
	void SetId(long id) { fId = id;}
	
	/** @brief Number of points to graphical resolution */
	virtual int NPoints(TPZGraphNode *n) = 0;
	
	virtual int NElements() = 0;
	/** @brief Sets a ith graphical node */
	virtual void SetNode(long i,TPZGraphNode *n);
	
	/** @brief Set dx style for connectivity information */
	virtual void Connectivity(TPZDrawStyle st = EDXStyle) = 0;
	/** @brief Draw coordinates of the graphical node */
	void DrawCo(TPZGraphNode *n, TPZDrawStyle st);
	
	/** @brief Draw solution of the graphical node */
	void DrawSolution(TPZGraphNode *n,TPZBlock<REAL> &Sol, TPZDrawStyle st);
	void DrawSolution(TPZGraphNode *n,int solind, TPZDrawStyle st);
	void DrawSolution(TPZGraphNode *n,TPZVec<int> &solind, TPZDrawStyle st);
	
	/** @brief Print the information of the graphical element */
	void Print(std::ostream &out);
	
	/** @brief Number of equations */
	virtual long EqNum(TPZVec<int> &co) = 0;
	
	
protected:
	/** @brief Computational element associated with graphical element */
	TPZCompEl *fCompEl;
	/** @brief Graphical mesh associated with graphical element */
	TPZGraphMesh *fGraphMesh;
	
	virtual void FirstIJ(int connect, TPZVec<int> &co, int &incr) = 0;
	
	virtual void NextIJ(int connect, TPZVec<int> &co, int incr) = 0;
	
	/** @brief Returns the number of the graphical node in the vector of connects */
	int ConnectNum(TPZGraphNode *n);
	
protected:
	/** @brief Id of the graphical element */
	long fId;
	
	/** @brief This method maps the index of a point to parameter space as a function of the number of divisions */
	virtual void QsiEta(TPZVec<int> &i, int imax, TPZVec<REAL> &qsieta);
	
};

#endif
