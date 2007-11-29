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

virtual ~TPZNonLinearAnalysis(void);

/**
 * It process a Newton's method to solve the non-linear problem.
 * It has the possibility of line search with parameter linesearch = true.
 */
void IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch = false);

/** Implements a golden section line search. */
void LineSearch(TPZFMatrix &Wn, TPZFMatrix &DeltaW, TPZFMatrix &NextW, REAL tol, int niter);

REAL SolutionNorm();

void ComputeTangent(TPZFMatrix &tangent, TPZVec<REAL> &coefs, int icase);

int NumCases();

void Residual(TPZFMatrix &residual, int icase);

void LoadSolution(TPZFMatrix &state);

void LoadState(TPZFMatrix &state);

};

#endif

