//$Id: pzelmat.h,v 1.2 2003-11-05 16:02:21 tiago Exp $

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
  /**block structure associated with fMat*/
  TPZBlock		*fBlock;
  /**pointer to a blocked matrix object*/
  TPZFMatrix		*fMat;
  /**vector of all nodes connected to the element*/
  TPZStack<int>	fConstrConnect;
  /**block structure associated with fConstrMat*/
  TPZBlock		*fConstrBlock;
  /**pointer to the constrained matrix object*/
  TPZFMatrix		*fConstrMat;

  TPZElementMatrix(int size) : fConnect(), fConstrConnect()
    {
      // allocating space for the vector of pointers
      fMat = NULL; // no space allocation for the matrix object
      fBlock = NULL;
      fConstrBlock = NULL;
      fConstrMat = NULL;

    }

  TPZElementMatrix() : fConnect(), fConstrConnect()
    {
      // allocating space for the vector of pointers
      fMat = NULL; // no space allocation for the matrix object
      fBlock = NULL;
      fConstrBlock = NULL;
      fConstrMat = NULL;

    }

  ~TPZElementMatrix(){
    // deletion of the vector does not delete the nodes (fortunately)
    if(fMat) delete fMat; // deletion of the matrix object
    if(fBlock) delete fBlock;
    if(fConstrBlock) delete fConstrBlock;
    if(fConstrMat) delete fConstrMat;
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
