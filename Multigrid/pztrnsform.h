#ifndef TTRANSFORMH
#define TTRANSFORMH

#include "pzfmatrix.h"

template<class T>
class TPZVec;

class TPZTransform {

  int fRow,fCol;     //Matrix dimensions
  TPZFMatrix fMult;	// multiplication matrix
  TPZFMatrix fSum;		// matrix used to sum
  REAL fStore[12];	// storage the matrix objects use to avoid
  // dynamic memory allocation
 public:

  TPZTransform(int dim);//square matrix

  TPZTransform(int fRow,int fCol);

  TPZTransform(const TPZTransform &tr);

  ~TPZTransform();

  TPZTransform &operator=(const TPZTransform &t);

  TPZFMatrix  & Mult() {return fMult;}

  TPZFMatrix  & Sum() {return fSum;}

  /**Sets the transformation matrices*/
  void SetMatrix(TPZFMatrix &mult,TPZFMatrix &sum);

  /**Multiply the transformation object to the right with right*/
  TPZTransform Multiply(TPZTransform &right);

  /**Transforms the vector*/
  void Apply(TPZVec<REAL> &vectorin,TPZVec<REAL> &vectorout);

  void PrintInputForm(ostream &out);

  int Compare(TPZTransform &t,REAL tol = 1.e-6);

};

#endif
