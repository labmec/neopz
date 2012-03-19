#ifndef CONVHYPER
#define CONVHYPER

#include "pzmathyperelastic.h"
#include "pzfmatrix.h"
#include "pzvec.h"

/// class used to check the consistency of the hyperelastic material
class TPZConvHyper : public TPZMatHyperElastic {

  TPZFMatrix<REAL> fState;
  TPZFMatrix<REAL> fPhi;
  TPZFMatrix<REAL> fDphi;
  TPZFMatrix<REAL> fAxes;
  TPZVec<REAL> fSol;
  TPZVec<REAL> fX;
  int fNumNod;

public :

  TPZConvHyper(int numnod, int mat, REAL E, REAL mu, REAL nu = -1., REAL lambda = -1., REAL coef1 = -1., REAL coef2 = -1., REAL coef3 = -1.);

  void ComputeTangent(TPZFMatrix<REAL> &tangent, TPZVec<REAL> &coefs, int icase);

  int NumCases();

  void Residual(TPZFMatrix<REAL> &residual, int icase);

  void LoadSolution(TPZFMatrix<REAL> &state);

};

#endif
