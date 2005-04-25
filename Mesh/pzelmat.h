//$Id: pzelmat.h,v 1.5 2005-04-25 02:31:48 phil Exp $

#ifndef ELMATHPP
#define ELMATHPP


#include "pzmatrix.h"
#include "pzblock.h"
#include "pzconnect.h"
#include "pzfmatrix.h"
#include "pzmanvector.h"
#include "pzstack.h"


/// this class associates an element matrix with the coeficients of its contribution in the global stiffness matrix
/**
This class groups all information associated with an element stiffness matrix so that it can be used independent of the element object itself
Objects of this class provide storage as well for the constrained stiffness matrix, i.e. the stiffness matrix from which the constrained connects have been eliminated
In future versions, the computation of the contraints will be incorporated in a method of this class
@ingroup interpolation
*/
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

  void Print(TPZCompMesh &mesh, std::ostream &out);

  void SetMatrixSize(short NumBli, short NumBlj, short BlSizei, short BlSizej);

  void SetMatrixMinSize(short NumBli, short NumBlj, short BlMinSizei, short BlMinSizej);

};

#endif
