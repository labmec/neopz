#ifndef MVGRIDH
#define MVGRIDH

#include "pzgraphmesh.h"
#include "pzvec.h"


class TPZBlock;

class TPZMVGraphMesh : public TPZGraphMesh {

public:

TPZMVGraphMesh(TPZCompMesh *cmesh, int dimension, TPZMaterial *mat);
TPZMVGraphMesh(TPZCompMesh *cmesh,int dim,TPZMVGraphMesh *graph,TPZMaterial *mat=0);

virtual void DrawMesh(int numcases);

virtual void DrawNodes();
virtual void DrawConnectivity();
virtual void DrawSolution(int step, REAL time,
						TPZVec<char *> &scalarnames, TPZVec<char *> &vectornames);
virtual void DrawSolution(TPZBlock &Sol);
virtual void DrawSolution(char *var = 0);

protected:
virtual void SequenceNodes();
	int fNumCases;
   int fNumSteps;

};

#endif

