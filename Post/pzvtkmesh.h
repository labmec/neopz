/**
 * @file
 * @brief Contains the TPZVTKGraphMesh class which implements the graphical mesh to VTK environment.
 */

#ifndef PZVTKMESH
#define PZVTKMESH

#include "pzgraphmesh.h"
#include "pzvec.h"

class TPZBlock;

/**
 * @ingroup post
 * @brief To export a graphical mesh to VTK environment. \ref post "Post processing"
 */
class TPZVTKGraphMesh : public TPZGraphMesh {
	
public:
	
    /** @brief Constructor for graphical mesh using VTK format with tensor variables */
    TPZVTKGraphMesh(TPZCompMesh *cmesh, int dimension, const std::set<int> & matids, const TPZVec<std::string> &scalnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensnames);
	/** @brief Copy constructor for graphical mesh using VTK format */
	TPZVTKGraphMesh(TPZCompMesh *cmesh,int dim,TPZVTKGraphMesh *graph);
	
	virtual void DrawMesh(int numcases) override;
	virtual void DrawNodes() override;
	virtual void DrawConnectivity(MElementType type) override;
	virtual void DrawSolution(int step, REAL time) override;
	virtual void DrawSolution(TPZBlock &Sol);
	virtual void DrawSolution(char *var = 0);
	
protected:
	virtual void SequenceNodes() override;
	int fNumCases;
	int fNumSteps;
	
};

#endif
