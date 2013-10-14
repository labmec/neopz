/**
 * @file
 * @brief Contains the TPZGraphNode class which implements the graphical node.
 */

#ifndef GRAFNODEH
#define GRAFNODEH

#include "pzconnect.h"
#include "pzgraphmesh.h"
#include "pzvec.h"

#include <iostream>

class TPZGraphMesh;
class TPZGraphEl;
template<class TVar>
class TPZBlock;

/**
 * @ingroup post
 * @brief To export a graphical node. \ref post "Post processing"
 */
class TPZGraphNode {
	
public:
	/** @brief Default constructor */
	TPZGraphNode();
	/** @brief Constructor for graphical node */
	TPZGraphNode(TPZConnect *cn, TPZGraphMesh *gm);
	/** @brief Simple destructor */
	~TPZGraphNode(void);
	
	//int ElIndex();
	long SequenceNumber() {return fSequenceNumber;}
	void SetSequenceNumber(long seqnum) {fSequenceNumber = seqnum;}
	void SetElement(TPZGraphEl *gel);
	void SetConnect(TPZConnect *connect);
	void SetGraphMesh(TPZGraphMesh *mesh);
	int NPoints();
	void SetPointNumber(long num);
	/** @brief Draw coordinates of the graphical node */
	void DrawCo(TPZDrawStyle st = EDXStyle);
	/** @brief Draw solution on the current connect for solutionid variable */
	void DrawSolution(int solutionid, TPZDrawStyle st = EDXStyle);
	void DrawSolution(TPZVec<int> &solutionid, TPZDrawStyle st= EDXStyle);
	void DrawSolution(TPZBlock<REAL> &sol, TPZDrawStyle st = EDXStyle);
	
	long FirstPoint();
	
	void Print(std::ostream &out);
	
protected:
	/** @brief Connect associated with current graphical node */
	TPZConnect *fConnect;
	/** @brief Graphical mesh related */
	TPZGraphMesh *fGraphMesh;
	/** @brief Graphical element related */
	TPZGraphEl *fGraphEl;
	long fPointNum;
	
private:
	long fSequenceNumber;
};

#endif
