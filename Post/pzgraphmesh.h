
#ifndef GRAFGRIDH
#define GRAFGRIDH


#include <iostream>
using namespace std;

#include <string.h>
#include "pzcompel.h"
#include "pzadmchunk.h"
#include "pzvec.h"

class TPZGraphNode;
class TPZCompMesh;
class TPZGraphEl;
class TPZFMatrix;
class TPZBlock;

enum TPZDrawStyle {EDXStyle,EMVStyle,EV3DStyle};

class TPZGraphMesh{
public:
  TPZGraphMesh(TPZCompMesh *cm, int dimension, TPZMaterial *mat);
  virtual ~TPZGraphMesh(void);
	
  TPZGraphNode &FindNode(int side);
  TPZGraphEl *FindElement(int sid);
  TPZAdmChunkVector<TPZGraphEl *> &ElementList();//MElementType type
  TPZAdmChunkVector<TPZGraphNode> &NodeMap();
  long NPoints();
  long NElements(MElementType type);
  int Res() {return fResolution;}
  void SetMaterial(TPZMaterial *mat) {fMaterial = mat;}
  virtual void SetCompMesh(TPZCompMesh *mesh, TPZMaterial *mat);
  ostream *Out();
  virtual void DrawNodes();
  virtual void DrawMesh(int numcases);
  virtual void DrawConnectivity(MElementType type);
  virtual void DrawSolution(TPZBlock &Sol);
  virtual void DrawSolution(char * var = 0);
  virtual void DrawSolution(int step, double time);
  TPZDrawStyle Style();
  void SetOutFile(ostream &out);
  void SetResolution(int res){ fResolution = res; SequenceNodes();}
	
  void Print(ostream &out);
  void SetNames(TPZVec<char *>&scalarnames, TPZVec<char *>&vecnames);
	
protected:
  TPZCompMesh *fCompMesh;
  TPZMaterial *fMaterial;
  int fDimension;
  TPZAdmChunkVector<TPZGraphEl *> fElementList;
  TPZAdmChunkVector<TPZGraphNode> fNodeMap;
  int fResolution;
  TPZDrawStyle fStyle;
  ostream *fOutFile;
  TPZVec<char *> fScalarNames, fVecNames;
  virtual void SequenceNodes();

  TPZInterpolatedElement *FindFirstInterpolatedElement(TPZCompMesh *mesh,int dimension);
	
public:
  virtual TPZMaterial *Material();
  virtual TPZCompMesh *Mesh() { return fCompMesh;}
};

inline TPZDrawStyle TPZGraphMesh::Style(){
  return fStyle;
}


#endif

