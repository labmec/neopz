//$Id: pzelmat.h,v 1.3 2004-04-05 14:09:35 phil Exp $

#ifndef ELMATHPP
#define ELMATHPP


#include "pzmatrix.h"
#include "pzblock.h"
#include "pzconnect.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzstack.h"


struct TPZElementMatrix {

  /**vector of pointers to TPZConnect objects*/
  TPZStack<int>	fConnect;
  /**pointer to a blocked matrix object*/
  TPZFNMatrix<1000>	fMat;
  /**block structure associated with fMat*/
  TPZBlock	       fBlock;
  /**vector of all nodes connected to the element*/
  TPZStack<int>	fConstrConnect;
  /**pointer to the constrained matrix object*/
  TPZFNMatrix<1000>		fConstrMat;
  /**block structure associated with fConstrMat*/
  TPZBlock		fConstrBlock;

  TPZElementMatrix() : fConnect(), fMat(0,0), fBlock(&fMat),  fConstrConnect(), fConstrMat(0,0), fConstrBlock(&fConstrMat)
    {
    }


  ~TPZElementMatrix(){
  }

  /**returns the number of nodes of TElementMatrix*/
  int NConnects(){
    return fConnect.NElements();
  }

  /**returns the pointer to the ith node of the element*/
  int ConnectIndex(int i){
    return fConnect[i];
  }

  void Print(TPZCompMesh &mesh, ostream &out);

  void SetMatrixSize(short NumBli, short NumBlj, short BlSizei, short BlSizej);

  void SetMatrixMinSize(short NumBli, short NumBlj, short BlMinSizei, short BlMinSizej);

};

#endif
