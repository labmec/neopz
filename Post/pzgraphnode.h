/**
 * @file
 * @brief Contains the TPZGraphNode class which implements the graphical node.
 */

#ifndef GRAFNODEH
#define GRAFNODEH

#include "pzconnect.h"
#include "pzvec.h"
#include "TPZDrawStyle.h"

#include <iostream>

class TPZGraphMesh;
class TPZGraphEl;
class TPZBlock;

/**
 * @ingroup post
 * @brief To export a graphical node. \ref post "Post processing"
 */
class TPZGraphNode : public TPZSavable {
	
public:
	/** @brief Default constructor */
	TPZGraphNode();
	/** @brief Constructor for graphical node */
	TPZGraphNode(TPZConnect *cn, TPZGraphMesh *gm);
	/** @brief Simple destructor */
	~TPZGraphNode(void);
            int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;

	//int ElIndex();
	int64_t SequenceNumber() {return fSequenceNumber;}
	void SetSequenceNumber(int64_t seqnum) {fSequenceNumber = seqnum;}
	void SetElement(TPZGraphEl *gel);
	void SetConnect(TPZConnect *connect);
	void SetGraphMesh(TPZGraphMesh *mesh);
	int NPoints();
	void SetPointNumber(int64_t num);
	/** @brief Draw coordinates of the graphical node */
	void DrawCo(TPZDrawStyle st = EDXStyle);
	/** @brief Draw solution on the current connect for solutionid variable */
	void DrawSolution(int solutionid, TPZDrawStyle st = EDXStyle);
	void DrawSolution(TPZVec<int> &solutionid, TPZDrawStyle st= EDXStyle);
	void DrawSolution(TPZBlock &sol, TPZDrawStyle st = EDXStyle);
	
	int64_t FirstPoint();
	
	void Print(std::ostream &out);
	
protected:
	/** @brief Connect associated with current graphical node */
	TPZConnect *fConnect;
	/** @brief Graphical mesh related */
	TPZGraphMesh *fGraphMesh;
	/** @brief Graphical element related */
	TPZGraphEl *fGraphEl;
	int64_t fPointNum;
	
private:
	int64_t fSequenceNumber;
};

#endif
