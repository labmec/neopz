#ifndef DXMESHH
#define DXMESHH

#include "pzgraphmesh.h"
#include "pzstack.h"


class TPZDXGraphMesh : public TPZGraphMesh {

  int fNextDataField;
  int fNumConnectObjects[8];
  int fElConnectivityObject[8];
  int fNodePosObject[8];
  int fNormalObject;
  TPZStack<REAL> fTimes;
  TPZStack<int> fFirstFieldValues[3];
  int fNumCases;
  char *fElementType;

public:

  TPZDXGraphMesh(TPZCompMesh *mesh, int dimension, TPZMaterial *mat, TPZVec<char *> &scalarnames, TPZVec<char *> &vecnames);
  TPZDXGraphMesh(TPZCompMesh *cmesh,int dim,TPZDXGraphMesh *graph,TPZMaterial *mat=0);

  virtual ~TPZDXGraphMesh();


  virtual void DrawMesh(int numcases);
  virtual void DrawSolution(TPZBlock &Sol);
  virtual void DrawSolution(char * var = 0);
  virtual void DrawSolution(int step, REAL time);

  int  NNodes();
  char *ElementName();
  void Close();
	
	/** Jorge 16/06/2001 */
	/** To redefinition of the graph mesh */
  void SetNumCases(int numcases) { fNumCases = numcases; }

private:
  void DrawNormals(int numnormals);
};

#endif
