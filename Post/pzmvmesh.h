#ifndef MVGRIDH
#define MVGRIDH

#include "pzgraphmesh.h"
#include "pzvec.h"


class TPZBlock;

/**
 * @ingroup post
 * @brief Implements graphical mesh to MVGraph package.
 */
/** 
 * MVGraph: Multivariate Interactive Visualization can to be obtained from: \n
 * \t http://cran.r-project.org/web/packages/mvgraph/index.html
 */
class TPZMVGraphMesh : public TPZGraphMesh {

public:

	TPZMVGraphMesh(TPZCompMesh *cmesh, int dimension, TPZAutoPointer<TPZMaterial> mat);
TPZMVGraphMesh(TPZCompMesh *cmesh,int dim,TPZMVGraphMesh *graph,TPZAutoPointer<TPZMaterial> mat);

virtual void DrawMesh(int numcases);

virtual void DrawNodes();
virtual void DrawConnectivity(MElementType type);
virtual void DrawSolution(int step, REAL time);
virtual void DrawSolution(TPZBlock &Sol);
virtual void DrawSolution(char *var = 0);

protected:
virtual void SequenceNodes();
	int fNumCases;
   int fNumSteps;

};

#endif

