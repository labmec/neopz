#ifndef V3DGRAFGH
#define V3DGRAFGH

#include "pzgraphmesh.h"
#include "pzvec.h"

//extern template class TPZVec<char *>;

class TPZV3DGraphMesh : public TPZGraphMesh {

public:

	TPZV3DGraphMesh(TPZCompMesh *cmesh, int dimension, TPZMaterial *mat);
  TPZV3DGraphMesh(TPZCompMesh *cmesh,int dim,TPZV3DGraphMesh *graph,TPZMaterial *mat = 0);

virtual void DrawMesh(int numcases);
	// Draw the nodal coordinates and the connectivity

virtual void DrawSolution(TPZBlock &Sol);
	// Draw the solution associated with Sol (not implemented)
virtual void DrawSolution(char * var = 0);
	// Draw the solution associated with the variable name
virtual void DrawSolution(int step, double time,
						TPZVec<char *> &scalarnames, TPZVec<char *> &vectornames);
	// Draw the solution sequence
virtual void SequenceNodes();

private :
	TPZCompMesh *fMesh;
	int fNumCases;
	int fInterval;
   int fLoadStep;
	int fTotScal;
	int fNumScal[6];
};

#endif

