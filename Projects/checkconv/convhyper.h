#ifndef CONVHYPER
#define CONVHYPER

#include "pzmathyperelastic.h"
#include "pzfmatrix.h"
#include "pzvec.h"

/// class used to check the consistency of the hyperelastic material
class TPZConvHyper : public TPZMatHyperElastic {

  TPZFMatrix fState;
  TPZFMatrix fPhi;
  TPZFMatrix fDphi;
  TPZFMatrix fAxes;
  TPZVec<REAL> fSol;
  TPZVec<REAL> fX;
  int fNumNod;

public :

  TPZConvHyper(int numnod, int mat, REAL E, REAL mu, REAL nu = -1., REAL lambda = -1., REAL coef1 = -1., REAL coef2 = -1., REAL coef3 = -1.);

  void ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &coefs, int icase);

  int NumCases();

  void Residual(TPZFMatrix &residual, int icase);

  void LoadSolution(TPZFMatrix &state);

};

#endif
