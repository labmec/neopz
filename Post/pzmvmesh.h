/**
 * @file
 * @brief Contains the TPZMVGraphMesh class which implements graphical mesh to MVGraph package.
 */

#ifndef MVGRIDH
#define MVGRIDH

#include "pzgraphmesh.h"
#include "pzvec.h"

class TPZBlock;

/**
 * @ingroup post
 * @brief Implements graphical mesh to MVGraph package. \ref post "Post processing"
 */
/** 
 * MVGraph: Multivariate Interactive Visualization can to be obtained from: \n
 * <a href="http://cran.r-project.org/web/packages/mvgraph/index.html">Lattes</a>
 */
class TPZMVGraphMesh : public TPZGraphMesh {
	
public:
	
	/** @brief Constructor for graphical mesh using MVGraph format */
    TPZMVGraphMesh(TPZCompMesh *cmesh, int dimension, const std::set<int> & matids, const TPZVec<std::string> &scalarnames, const TPZVec<std::string> &vecnames);
	/** @brief Copy constructor for graphical mesh using MVGraph format */
	TPZMVGraphMesh(TPZCompMesh *cmesh,int dim,TPZMVGraphMesh *graph);
	
	/** @brief Draw graphical mesh */
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

