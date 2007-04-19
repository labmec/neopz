//$Id: pzelmat.h,v 1.8 2007-04-19 12:21:36 tiago Exp $

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

  enum MType{EF = 1, EK = 2};

  MType fType;

  TPZCompMesh * fMesh;

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
  
  TPZManVector<int> fDestinationIndex, fSourceIndex;
  
  int fNumStateVars;
  
  ///Reset the data structure
  void Reset()
  {
    fConnect.Resize(0);
    fMat.Resize(0,0);
    fBlock.SetNBlocks(0);
    fConstrConnect.Resize(0);
    fConstrMat.Resize(0,0);
    fConstrBlock.SetNBlocks(0);
  }

  TPZElementMatrix(TPZCompMesh *mesh, MType type) : fType(type), fMesh(mesh), fConnect(), fMat(0,0), fBlock(&fMat),  fConstrConnect(), fConstrMat(0,0), fConstrBlock(&fConstrMat), fNumStateVars(0)
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

  void Print(std::ostream &out);

  void SetMatrixSize(short NumBli, short NumBlj, short BlSizei, short BlSizej);

  void SetMatrixMinSize(short NumBli, short NumBlj, short BlMinSizei, short BlMinSizej);

  void ComputeDestinationIndices();


  /**
   * Returns true if the element has at least one dependent node
   * returns false otherwise
   */
  bool HasDependency();

  /**
   * Apply the constraints applied to the nodes by transforming the tangent
   * matrix and right hand side
   * @param ekmat element stiffness matrix
   * @param efmat element loads matrix
   */
  void ApplyConstraints();

};

#endif
