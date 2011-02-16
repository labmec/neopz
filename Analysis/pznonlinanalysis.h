#ifndef NONLINANALYSISH
#define NONLINANALYSISH

#include "pzanalysis.h"
#include "pzcompel.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include <iostream>
class TPZCompMesh;
class TPZFMatrix;

class TPZNonLinearAnalysis : public TPZAnalysis {

public:

TPZNonLinearAnalysis(TPZCompMesh *mesh,std::ostream &out);

TPZNonLinearAnalysis();

virtual ~TPZNonLinearAnalysis();

/**
 * It process a Newton's method to solve the non-linear problem.
 * It has the possibility of line search with parameter linesearch = true.
 */
virtual void IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch = false, bool checkconv = false);

/** Implements a golden section line search. */
void LineSearch(const TPZFMatrix &Wn, TPZFMatrix DeltaW, TPZFMatrix &NextW, REAL tol, int niter);

REAL SolutionNorm();

virtual void ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &coefs, int icase);

int NumCases();

virtual void Residual(TPZFMatrix &residual, int icase);

virtual void LoadSolution(const TPZFMatrix &state);

void LoadState(TPZFMatrix &state);

};

#endif

