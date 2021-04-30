/**
 * @file
 * @brief Contains the TPZGraphEl1d class which implements the graphical one dimensional element.
 */

#ifndef GRAFEL1DH
#define GRAFEL1DH

#include "pzgraphel.h"
#include "pzvec.h"

class TPZGraphMesh;
class TPZGraphNode;
class TPZBlock;

/**
 * @ingroup post
 * @brief To export a graphical one dimensional element. \ref post "Post processing"
 */
class TPZGraphEl1d : public TPZGraphEl
{
public:
	/** @brief Constructor for graphical element to computational one dimensional element */
	TPZGraphEl1d(TPZCompEl *ce, TPZGraphMesh *gg);
	
	virtual int NConnects() override { return 3;}
	
	virtual MElementType Type() override {return EOned;}
	
	virtual int ExportType(TPZDrawStyle st) override;
	
	virtual int NNodes() override;
	
	virtual TPZGraphNode *Connect(int64_t i) override {return fConnects[i];}
	
	virtual int NPoints(TPZGraphNode *n) override;
	
	virtual int NElements() override;
	
    /** @brief the parametric dimension of the element */
    virtual int Dimension() override
    {
        return 1;
    }

    virtual void Connectivity(TPZDrawStyle st = EDXStyle) override;
	
	void Print(std::ostream &out);
	
	virtual int64_t EqNum( TPZVec<int> &co) override;
	
	
protected:
	
	virtual void FirstIJ(int no, TPZVec<int> &co, int &incr) override;
	
	virtual void NextIJ(int no,  TPZVec<int> &co, int incr) override;
	
	/** @brief Graphical nodes vector (by connect of the computational element) */
	TPZGraphNode *fConnects[3];
	
	virtual void SetNode(int64_t i,TPZGraphNode *gno) override {
		fConnects[i] = gno;
	}
	
};

#endif
