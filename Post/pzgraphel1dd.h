/**
 * @file
 * @brief Contains the TPZGraphEl1dd class which implements the graphical one dimensional discontinuous element.
 */

#ifndef GRAFEL1DDH
#define GRAFEL1DDH


#include "pzgraphel.h"
#include "pzvec.h"

class TPZGraphMesh;
class TPZGraphNode;
class TPZBlock;

/**
 * @ingroup post
 * @brief To export a graphical one dimensional discontinuous element. \ref post "Post processing"
 */
class TPZGraphEl1dd : public TPZGraphEl
{
public:
	/** @brief Constructor for graphical element to computational one dimensional discontinuous element */
	TPZGraphEl1dd(TPZCompEl *ce, TPZGraphMesh *gg);
	
	virtual int NPoints(TPZGraphNode *n) override;
	
	virtual int NElements() override;
	
	virtual void Connectivity(TPZDrawStyle st = EDXStyle) override;
	
	void Print(std::ostream &out);
	
	virtual int64_t EqNum(TPZVec<int> &co) override;
	
	virtual int ExportType(TPZDrawStyle st) override;
	
	virtual int NNodes() override;
	
    /** @brief the parametric dimension of the element */
    virtual int Dimension() override
    {
        return 1;
    }

	
protected:
	
	virtual void FirstIJ(int no, TPZVec<int> &co, int &incr) override;
	
	virtual void NextIJ(int no, TPZVec<int> &co, int incr) override;
	
	/** @brief Graphical node (connect) to discontinuous graphical element */
	TPZGraphNode *fConnect;
	
	virtual void SetNode(int64_t i,TPZGraphNode *gno) override {
		fConnect = gno;
	}
	
	int NConnects() override { return 1; }
	MElementType Type() override {return EOned;}
	
	TPZGraphNode *Connect(int64_t i) override {return fConnect;}
	
};

#endif
